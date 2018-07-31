
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#ifndef HEMELB_NET_MPICOMMUNICATOR_HPP
#define HEMELB_NET_MPICOMMUNICATOR_HPP

#include "net/MpiDataType.h"
#include "net/MpiConstness.h"

namespace hemelb
{
  namespace net
  {
    template<typename T>
    void MpiCommunicator::Broadcast(T& val, const int root) const
    {
      HEMELB_MPI_CALL(
          MPI_Bcast,
          (&val, 1, MpiDataType<T>(), root, *this)
      );
    }
    template<typename T>
    void MpiCommunicator::Broadcast(std::vector<T>& vals, const int root) const
    {
      HEMELB_MPI_CALL(
          MPI_Bcast,
          (&vals[0], vals.size(), MpiDataType<T>(), root, *this)
      );
    }

    template<typename T>
    T MpiCommunicator::AllReduce(const T& val, const MPI_Op& op) const
    {
      T ans;
      HEMELB_MPI_CALL(
          MPI_Allreduce,
          (MpiConstCast(&val), &ans, 1, MpiDataType<T>(), op, *this)
      );
      return ans;
    }

    template<typename T>
    std::vector<T> MpiCommunicator::AllReduce(const std::vector<T>& vals, const MPI_Op& op) const
    {
      std::vector<T> ans(vals.size());
      HEMELB_MPI_CALL(
          MPI_Allreduce,
          (MpiConstCast(&vals[0]), &ans[0], vals.size(), MpiDataType<T>(), op, *this)
      );
      return ans;
    }

    template<typename T>
    T MpiCommunicator::Reduce(const T& val, const MPI_Op& op, const int root) const
    {
      T ans;
      HEMELB_MPI_CALL(
          MPI_Reduce,
          (MpiConstCast(&val), &ans, 1, MpiDataType<T>(), op, root, *this)
      );
      return ans;
    }

    template<typename T>
    std::vector<T> MpiCommunicator::Reduce(const std::vector<T>& vals, const MPI_Op& op,
                                           const int root) const
    {
      std::vector<T> ans;
      T* recvbuf = NULL;

      if (Rank() == root)
      {
        // Standard says the address of receive buffer only matters at the root.
        ans.resize(vals.size());
        recvbuf = &ans[0];
      }

      HEMELB_MPI_CALL(
          MPI_Reduce,
          (MpiConstCast(&vals[0]), recvbuf, vals.size(), MpiDataType<T>(), op, root, *this)
      );
      return ans;
    }

    template<typename T>
    std::vector<T> MpiCommunicator::Gather(const T& val, const int root) const
    {
      std::vector<T> ans;
      T* recvbuf = NULL;

      if (Rank() == root)
      {
        // Standard says the address of receive buffer only matters at the root.
        ans.resize(Size());
        recvbuf = &ans[0];
      }
      HEMELB_MPI_CALL(
          MPI_Gather,
          (MpiConstCast(&val), 1, MpiDataType<T>(),
              recvbuf, 1, MpiDataType<T>(),
              root, *this)
      );
      return ans;
    }

    template <typename T>
    T MpiCommunicator::Scatter(const std::vector<T>& vals, const int root) const {
      T ans;
      const T* ptr = (Rank() == root) ? vals.data() : nullptr;
      HEMELB_MPI_CALL(
		      MPI_Scatter,
		      (MpiConstCast(ptr), 1, MpiDataType<T>(),
		       &ans, 1, MpiDataType<T>(),
		       root, *this)
		      );
      return ans;
    }

    template <typename T>
    std::vector<T> MpiCommunicator::Scatter(const std::vector<T>& vals, const size_t n, const int root) const {
      std::vector<T> ans(n);
      const T* ptr = (Rank() == root) ? vals.data() : nullptr;
      HEMELB_MPI_CALL(
		      MPI_Scatter,
		      (MpiConstCast(ptr), n, MpiDataType<T>(),
		       ans.data(), n, MpiDataType<T>(),
		       root, *this)
		      );
      return ans;
    }

    template<typename T>
    std::vector<T> MpiCommunicator::AllGather(const T& val) const
    {
      std::vector<T> ans(Size());
      T* recvbuf =  &ans[0];

      HEMELB_MPI_CALL(
          MPI_Allgather,
          (MpiConstCast(&val), 1, MpiDataType<T>(),
              recvbuf, 1, MpiDataType<T>(),
              *this)
          );
      return ans;
    }

    template <typename T>
    std::vector<T> MpiCommunicator::AllToAll(const std::vector<T>& vals) const
    {
      std::vector<T> ans(vals.size());
      HEMELB_MPI_CALL(
          MPI_Alltoall,
          (MpiConstCast(&vals[0]), 1, MpiDataType<T>(),
           &ans[0], 1, MpiDataType<T>(),
           *this)
      );
      return ans;
    }

    template <typename T>
    void MpiCommunicator::Send(const T& val, int dest, int tag) const
    {
      HEMELB_MPI_CALL(
          MPI_Send,
          (MpiConstCast(&val), 1, MpiDataType<T>(), dest, tag, *this)
      );
    }
    template <typename T>
    void MpiCommunicator::Send(const std::vector<T>& vals, int dest, int tag) const
    {
      HEMELB_MPI_CALL(
          MPI_Send,
          (MpiConstCast(&vals[0]), vals.size(), MpiDataType<T>(), dest, tag, *this)
      );
    }

    template <typename T>
    void MpiCommunicator::Receive(T& val, int src, int tag, MPI_Status* stat) const
    {
      HEMELB_MPI_CALL(
          MPI_Recv,
          (&val, 1, MpiDataType<T>(), src, tag, *this, stat)
      );
    }
    template <typename T>
    void MpiCommunicator::Receive(std::vector<T>& vals, int src, int tag, MPI_Status* stat) const
    {
      HEMELB_MPI_CALL(
          MPI_Recv,
          (&vals, vals.size(), MpiDataType<T>(), src, tag, *this, stat)
      );
    }
  }
}

#endif // HEMELB_NET_MPICOMMUNICATOR_HPP
