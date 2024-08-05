// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_NET_IOCOMMUNICATOR_H
#define HEMELB_NET_IOCOMMUNICATOR_H

#include "net/MpiCommunicator.h"

namespace hemelb::net {
    /**
     * An MPI communicator, which has two special sub communicators:
     *
     * - one for each MPI shared mem region (i.e. per node)
     *
     * - one made up of all the node-rank-zero processes (i.e. size of
     * number of nodes, one member per node)
     *
     * The top level communicator also has a special IO rank.
     */
    class IOCommunicator : public MpiCommunicator
    {
        MpiCommunicator nodeComm;
        bool amNodeLeader;
        MpiCommunicator leadersComm;
    public:
        static constexpr int IO_RANK = 0;

        IOCommunicator() = default;
        explicit IOCommunicator(const MpiCommunicator& comm);

        inline bool OnIORank() const {
            return Rank() == IO_RANK;
        }
        constexpr int GetIORank() const {
            return IO_RANK;
        }

        inline MpiCommunicator const& GetNodeComm() const {
            return nodeComm;
        }
        inline bool AmNodeLeader() const {
            return amNodeLeader;
        }
        inline MpiCommunicator const& GetLeadersComm() const {
            return leadersComm;
        }
    };
}

#endif /* HEMELB_NET_IOCOMMUNICATOR_H */
