// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_NET_SPARSEEXCHANGE_H
#define HEMELB_NET_SPARSEEXCHANGE_H

#include <span>
#include <utility>
#include "net/MpiCommunicator.h"

namespace hemelb::net
{
    // Do a *single* collective sparse all to all of possibly variable size.
    //
    // The data swapped has to be of a single type `T`.
    // For variable message size use std::dynamic_extent (like span).
    //
    // Each process can make as many sends to whichever processes it like.
    //
    // Then all processes must collectively call receive with appropriate
    // callbacks (see below).
    //
    // The object is then only fit to be destroyed.
    template <typename T, std::size_t MSG_SIZE = std::dynamic_extent>
    struct sparse_exchange {
        MpiCommunicator comm;
        int tag;

        MPI_Datatype dtype = MPI_DATATYPE_NULL;
        std::vector<MPI_Request> send_reqs;
        MpiRequest barrier_req;
        int my_sends_received = 0;
        bool all_sends_received = false;

        // Specify communicator and a tag to use.
        sparse_exchange(MpiCommunicator c, int t) :
                comm(std::move(c)), tag(t),
                dtype(MpiDataType<T>()) {
        }

        // Send a single value, if allowed
        void send(T const &val, int to_rank) requires (MSG_SIZE == 1 || MSG_SIZE == std::dynamic_extent) {
            MPI_Request r;
            MpiCall{MPI_Issend}(&val, 1, dtype, to_rank, tag, comm, &r);
            send_reqs.push_back(r);
        }

        // Send a range of values
        template <std::size_t N>
        requires (N == std::dynamic_extent || MSG_SIZE == std::dynamic_extent || N == MSG_SIZE)
        void send(std::span<T const, N> data, int to_rank) {
            // If the exchange's size is dynamic, can't check
            if constexpr (MSG_SIZE != std::dynamic_extent) {
                if constexpr (N == std::dynamic_extent) {
                    if (data.size() != MSG_SIZE)
                        throw (Exception() << "Wrong message size!");
                } // The else has been checked by the requires clause above!
            }

            MPI_Request r;
            MpiCall{MPI_Issend}(data.data(), std::ssize(data), dtype, to_rank, tag, comm, &r);
            send_reqs.push_back(r);
        }

        // Do the collective receive of messages.
        //
        // First arg decides what to do about the size of received messages.
        //
        // It is passed the source rank and the number of data items.
        //
        // It must return a pointer to a contiguous piece of memory which
        // can hold that number of items.
        //
        // The pointer must remain valid until the second argument is called
        // with the source rank and the pointer.
        template<std::invocable<int, int> SizeHandler, std::invocable<int, T*> RecvHandler>
        void receive(SizeHandler&& sizeHandler, RecvHandler&& recvHandler) {
            while (!all_sends_received) {
                // Check for arriving messages
                MPI_Status status;
                int msg_available;
                MpiCall{MPI_Iprobe}(MPI_ANY_SOURCE, tag, comm, &msg_available, &status);
                if (msg_available) {
                    // We have one - check it's correct and push data onto vec
                    int nrecv;
                    MpiCall{MPI_Get_count}(&status, dtype, &nrecv);
                    if constexpr (MSG_SIZE != std::dynamic_extent) {
                        if (nrecv != MSG_SIZE)
                            throw (Exception() << "Wrong message size!");
                    }
                    T *buf = sizeHandler(status.MPI_SOURCE, nrecv);
                    MpiCall{MPI_Recv}(buf, nrecv, dtype, status.MPI_SOURCE, tag, comm, MPI_STATUS_IGNORE);
                    recvHandler(status.MPI_SOURCE, buf);
                }

                if (!my_sends_received) {
                    // check to see if all this rank's sends have been received
		  MpiCall{MPI_Testall}(std::ssize(send_reqs), send_reqs.data(), &my_sends_received, MPI_STATUSES_IGNORE);
                    if (my_sends_received) {
                        // They have! Signal this to the rest of the communicator.
                        barrier_req = comm.Ibarrier();
                        // This Testall/Ibarrier won't run again now, but this rank
                        // will keep checking for incoming messages until all processes
                        // have started the barrier.
                    }
                } else {
                    all_sends_received = barrier_req.Test();
                }
            }
        }
    };
}

#endif
