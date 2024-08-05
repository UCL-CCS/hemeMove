// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#ifndef HEMELB_NET_MPICOMMUNICATOR_HPP
#define HEMELB_NET_MPICOMMUNICATOR_HPP

#include <numeric>
#include "net/MpiDataType.h"

namespace hemelb::net
{
    template<typename T>
    void MpiCommunicator::Broadcast(T& val, const int root) const
    {
      MpiCall{MPI_Bcast}(&val, 1, MpiDataType<T>(), root, *this);
    }
    template<typename T, std::size_t N>
    void MpiCommunicator::Broadcast(std::span<T, N> vals, const int root) const
    {
      MpiCall{MPI_Bcast}(vals.data(), vals.size(), MpiDataType<T>(), root, *this);
    }

    template<typename T>
    T MpiCommunicator::AllReduce(const T& val, const MPI_Op& op) const
    {
      T ans;
      HEMELB_MPI_CALL(MPI_Allreduce, (&val, &ans, 1, MpiDataType<T>(), op, *this));
      return ans;
    }

    template<typename T>
    void MpiCommunicator::AllReduceInPlace(const std::span<T>& vals, const MPI_Op& op) const
    {
        HEMELB_MPI_CALL(MPI_Allreduce,
                        (MPI_IN_PLACE, vals.data(), vals.size(), MpiDataType<T>(), op, *this));
    }

    template<typename T>
    std::vector<T> MpiCommunicator::AllReduce(const std::vector<T>& vals, const MPI_Op& op) const
    {
      std::vector<T> ans(vals.size());
      HEMELB_MPI_CALL(MPI_Allreduce,
                      (vals.data(), ans.data(), vals.size(), MpiDataType<T>(), op, *this));
      return ans;
    }

    template<typename T>
    T MpiCommunicator::Scan(const T& val, const MPI_Op& op) const
    {
      T ans;
      HEMELB_MPI_CALL(MPI_Scan, (&val, &ans, 1, MpiDataType<T>(), op, *this));
      return ans;
    }

    template<typename T>
    MpiRequest MpiCommunicator::Iexscan(T const& val, T& dest, const MPI_Op &op) const
    {
        MpiRequest ans;
        HEMELB_MPI_CALL(MPI_Iexscan, (&val, &dest, 1, MpiDataType<T>(), op, *this, &ans.req));
        return ans;
    }

    template<typename T>
    T MpiCommunicator::Reduce(const T& val, const MPI_Op& op, const int root) const
    {
      T ans;
      HEMELB_MPI_CALL(MPI_Reduce, (&val, &ans, 1, MpiDataType<T>(), op, root, *this));
      return ans;
    }
    template<typename T, std::size_t N>
    void MpiCommunicator::Reduce(std::span<T, N> dest, std::span<const T, N> vals, const MPI_Op& op, const int root) const
    {
        HEMELB_MPI_CALL(MPI_Reduce, (vals.data(), dest.data(), vals.size(), MpiDataType<T>(), op, root, *this));
    }

    template<typename T>
    std::vector<T> MpiCommunicator::Reduce(const std::vector<T>& vals, const MPI_Op& op,
                                           const int root) const
    {
      std::vector<T> ans;
      T* recvbuf = nullptr;

      if (Rank() == root)
      {
        // Standard says the address of receive buffer only matters at the root.
        ans.resize(vals.size());
        recvbuf = ans.data();
      }

      HEMELB_MPI_CALL(MPI_Reduce,
                      (vals.data(), recvbuf, vals.size(), MpiDataType<T>(), op, root, *this));
      return ans;
    }

    template<typename T>
    std::vector<T> MpiCommunicator::Gather(const T& val, const int root) const
    {
      std::vector<T> ans;
      T* recvbuf = nullptr;

      if (Rank() == root)
      {
        // Standard says the address of receive buffer only matters at the root.
        ans.resize(Size());
        recvbuf = ans.data();
      }
      HEMELB_MPI_CALL(MPI_Gather,
                      (&val, 1, MpiDataType<T>(), recvbuf, 1, MpiDataType<T>(), root, *this));
      return ans;
    }

    template<typename T>
    std::vector<T> MpiCommunicator::Gather(const std::vector<T>& val, const int root) const
    {
      std::vector<T> ans;
      std::vector<int> offsets;
      auto const counts = Gather(static_cast<int>(val.size()), root);

      if (Rank() == root)
      {
        offsets.push_back(0);
        for(auto const &c: counts)
        {
          offsets.push_back(offsets.back() + c);
        }

        auto const N = std::accumulate(counts.begin(), counts.end(), 0);
        ans.reserve(N+1);
        ans.resize(N);
      }

      HEMELB_MPI_CALL(MPI_Gatherv,
                      (val.data(), val.size(), MpiDataType<T>(),
                       Rank() == root ? ans.data(): nullptr,
                       Rank() == root ? counts.data(): nullptr,
                       Rank() == root ? offsets.data(): nullptr,
                       MpiDataType<T>(), root, *this));
      return ans;
    }

    template <typename T>
    T MpiCommunicator::Scatter(const std::vector<T>& vals, const int root) const {
      T ans;
      const T* ptr = (Rank() == root) ? vals.data() : nullptr;
      HEMELB_MPI_CALL(
		      MPI_Scatter,
		      (ptr, 1, MpiDataType<T>(),
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
		      (ptr, n, MpiDataType<T>(),
		       ans.data(), n, MpiDataType<T>(),
		       root, *this)
		      );
      return ans;
    }

    template<typename T>
    std::vector<T> MpiCommunicator::AllGather(const T& val) const
    {
      std::vector<T> ans(Size());

      HEMELB_MPI_CALL(MPI_Allgather,
                      (&val, 1, MpiDataType<T>(), ans.data(), 1, MpiDataType<T>(), *this));
      return ans;
    }

    template <typename T>
    displaced_data<T> MpiCommunicator::AllGatherV(const std::vector<T>& local_data) const {
        std::vector<int> per_rank_sizes = AllGather((int) local_data.size());
        auto ans = displaced_data<T>{per_rank_sizes};
        HEMELB_MPI_CALL(MPI_Allgatherv,
                        (local_data.data(), local_data.size(), MpiDataType<T>(),
                         ans.data.data(), per_rank_sizes.data(), ans.displacements.data(), MpiDataType<T>(),
                         *this)
        );
        return ans;
      }
    template<typename T>
    std::vector<T> MpiCommunicator::AllNeighGather(const T& val) const
    {
      std::vector<T> ans(GetNeighborsCount());

      HEMELB_MPI_CALL(MPI_Neighbor_allgather,
                      (&val, 1, MpiDataType<T>(), ans.data(), 1, MpiDataType<T>(), *this));

      return ans;
    }

    template<typename T>
    displaced_data<T> MpiCommunicator::AllNeighGatherV(const std::vector<T>& val) const
    {
      std::vector<int> valSizes = AllNeighGather((int) val.size());
      auto ans = displaced_data<T>{valSizes};
      HEMELB_MPI_CALL(
              MPI_Neighbor_allgatherv,
              (val.data(), val.size(), MpiDataType<T>(),
               ans.data.data(), valSizes.data(), ans.displacements.data(), MpiDataType<T>(),
               *this)
      );
      return ans;
    }

    template<typename T>
    std::vector<T> MpiCommunicator::AllToAll(const std::vector<T>& vals) const
    {
      std::vector<T> ans(vals.size());
      HEMELB_MPI_CALL(MPI_Alltoall,
                      (vals.data(), 1, MpiDataType<T>(), ans.data(), 1, MpiDataType<T>(), *this));
      return ans;
    }

    template<typename T>
    void MpiCommunicator::Send(const T& val, int dest, int tag) const
    {
      HEMELB_MPI_CALL(MPI_Send, (&val, 1, MpiDataType<T>(), dest, tag, *this));
    }
    template<typename T>
    void MpiCommunicator::Send(const std::vector<T>& vals, int dest, int tag) const
    {
      HEMELB_MPI_CALL(MPI_Send,
                      (vals.data(), vals.size(), MpiDataType<T>(), dest, tag, *this));
    }

    template <typename T>
    MpiRequest MpiCommunicator::Issend(std::span<T const> vals, int dest, int tag) const {
        MpiRequest ans;
        MpiCall{MPI_Issend}(vals.data(), vals.size(), MpiDataType<T>(), dest, tag, *commPtr, &ans.req);
        return ans;
    }
    template <typename T>
    MpiRequest MpiCommunicator::Issend(T const& val, int dest, int tag) const {
        MpiRequest ans;
        MpiCall{MPI_Issend}(&val, 1, MpiDataType<T>(), dest, tag, *commPtr, &ans.req);
        return ans;
    }

    template<typename T>
    void MpiCommunicator::Receive(T& val, int src, int tag, MPI_Status* stat) const
    {
      HEMELB_MPI_CALL(MPI_Recv, (&val, 1, MpiDataType<T>(), src, tag, *this, stat));
    }
    template<typename T>
    void MpiCommunicator::Receive(std::vector<T>& vals, int src, int tag, MPI_Status* stat) const
    {
      HEMELB_MPI_CALL(MPI_Recv, (vals.data(), vals.size(), MpiDataType<T>(), src, tag, *this, stat));
    }

    template <typename T>
    MpiRequest MpiCommunicator::Irecv(std::span<T> dest, int src, int tag) const {
        MpiRequest ans;
        MpiCall{MPI_Irecv}(dest.data(), dest.size(), MpiDataType<T>(), src, tag, *commPtr, &ans.req);
        return ans;
    }
    template <typename T>
    MpiRequest MpiCommunicator::Irecv(T &dest, int src, int tag) const {
        MpiRequest ans;
        MpiCall{MPI_Irecv}(&dest, 1, MpiDataType<T>(), src, tag, *commPtr, &ans.req);
        return ans;
    }
}

#endif // HEMELB_NET_MPICOMMUNICATOR_HPP
