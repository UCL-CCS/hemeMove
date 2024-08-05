// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <numeric>

#include "net/MpiCommunicator.h"
#include "net/MpiGroup.h"
#include "util/Iterator.h"

namespace hemelb::net
{
    void MpiRequest::Waitall(std::span<MpiRequest> reqs) {
        // These asserts check that there's no extra data in an
        // MpiRequest object, so we can effectively just pointer
        // alias a span of them.
        static_assert(sizeof(MpiRequest) == sizeof(MPI_Request));
        static_assert(alignof(MpiRequest)== alignof(MPI_Request));
        int N = std::ssize(reqs);
        MPI_Request* data = reqs.data() ? &reqs.data()->req : nullptr;
        MpiCall{MPI_Waitall}(N, data, MPI_STATUSES_IGNORE);
    }

    namespace {
      void Deleter(MPI_Comm* comm)
      {
        int finalized;
        HEMELB_MPI_CALL(MPI_Finalized, (&finalized));
        if (!finalized)
          HEMELB_MPI_CALL(MPI_Comm_free, (comm));
        delete comm;
      }
    }

    MpiCommunicator MpiCommunicator::World()
    {
      return {MPI_COMM_WORLD, false};
    }

    MpiCommunicator::MpiCommunicator() :
        commPtr(), communicatorSize(-1), localRankInCommunicator(-1)
    {
    }

    MpiCommunicator::MpiCommunicator(MPI_Comm communicator, bool owner) :
        commPtr(), communicatorSize(-1), localRankInCommunicator(-1)
    {
      if (communicator == MPI_COMM_NULL)
        return;

      if (owner)
      {
        commPtr.reset(new MPI_Comm(communicator), Deleter);
      }
      else
      {
        commPtr = std::make_shared<MPI_Comm>(communicator);
      }

      HEMELB_MPI_CALL(MPI_Comm_size, (communicator, &communicatorSize));
      HEMELB_MPI_CALL(MPI_Comm_rank, (communicator, &localRankInCommunicator));
    }

    MpiCommunicator::MpiCommunicator(mock_ctor_tag, int localRankInCommunicator, int communicatorSize) :
        commPtr(), communicatorSize(communicatorSize), localRankInCommunicator(localRankInCommunicator)
    {
    }

    bool operator==(const MpiCommunicator& comm1, const MpiCommunicator& comm2)
    {
      if (comm1)
      {
        if (comm2)
        {
          int result;
          HEMELB_MPI_CALL(MPI_Comm_compare, (comm1, comm2, &result));
          return result == MPI_IDENT;
        }
        return false;
      }
      return (!comm2);
    }

    bool operator!=(const MpiCommunicator& comm1, const MpiCommunicator& comm2)
    {
      return ! (comm1 == comm2);
    }

    void MpiCommunicator::Barrier() const
    {
      HEMELB_MPI_CALL(MPI_Barrier, (*commPtr));
    }

    void MpiCommunicator::Broadcast(std::string& val, const int root) const
    {
        auto len = val.size();
        Broadcast(len, root);
        if (Rank() != root) {
            val.resize(len);
        }
        Broadcast(std::span<char>{val.data(), len}, root);
    }

    MpiRequest MpiCommunicator::Ibarrier() const {
        MpiRequest ans;
        MpiCall{MPI_Ibarrier}(*commPtr, &ans.req);
        return ans;
    }

    MpiGroup MpiCommunicator::Group() const
    {
      MPI_Group grp;
      HEMELB_MPI_CALL(MPI_Comm_group, (*commPtr, &grp));
      return {grp, true};
    }

    MpiCommunicator MpiCommunicator::Create(const MpiGroup& grp) const
    {
      MPI_Comm newComm;
      HEMELB_MPI_CALL(MPI_Comm_create, (*commPtr, grp, &newComm));
      return {newComm, true};
    }

    void MpiCommunicator::Abort(int errCode) const
    {
      HEMELB_MPI_CALL(MPI_Abort, (*commPtr, errCode));
    }

    MpiCommunicator MpiCommunicator::Duplicate() const
    {
      MPI_Comm newComm;
      HEMELB_MPI_CALL(MPI_Comm_dup, (*commPtr, &newComm));
      return {newComm, true};
    }

    MpiCommunicator MpiCommunicator::DistGraphAdjacent(std::vector<int> my_neighbours, bool reorder) const
    {
        MPI_Comm newComm;
        HEMELB_MPI_CALL(MPI_Dist_graph_create_adjacent ,
                        (*commPtr,
                                my_neighbours.size(), my_neighbours.data(), MPI_UNWEIGHTED,
                                my_neighbours.size(), my_neighbours.data(), MPI_UNWEIGHTED,
                                MPI_INFO_NULL, reorder, &newComm)
        );
        return {newComm, true};
    }

    int MpiCommunicator::GetNeighborsCount() const
    {
        int n_in, n_out, weighted;
        HEMELB_MPI_CALL(MPI_Dist_graph_neighbors_count, (*commPtr, &n_in, &n_out, &weighted));
#ifndef NDEBUG
        if (weighted)
            throw (Exception() << "Only support unweighted graphs");
        if (n_in != n_out)
            throw (Exception() << "Only support bidirectional graphs");
#endif
        return n_in;
    }

    std::vector<int> MpiCommunicator::GetNeighbors() const
    {
        auto n = GetNeighborsCount();
        auto result = std::vector<int>(n);
        auto ignored = std::vector<int>(n);
        MPI_Dist_graph_neighbors(*commPtr, n, result.data(), MPI_UNWEIGHTED, n, ignored.data(), MPI_UNWEIGHTED);
#ifndef NDEBUG
        for (int i = 0; i < n; ++i)
            if (result[i] != ignored[i])
                throw (Exception() << "In and out connections differ");
#endif
        return result;
    }

    MpiCommunicator MpiCommunicator::Split(int color, int key) const
    {
      MPI_Comm newComm;
      HEMELB_MPI_CALL(MPI_Comm_split, (*commPtr, color, key, &newComm));
      return {newComm, true};
    }

    MpiCommunicator MpiCommunicator::SplitType(int type) const
    {
      MPI_Comm newComm;
      HEMELB_MPI_CALL(MPI_Comm_split_type, (*commPtr, type, Rank(), MPI_INFO_NULL, &newComm));
      return {newComm, true};
    }

    std::map<int, int> MpiCommunicator::RankMap(MpiCommunicator const &valueComm) const
    {
      MPI_Group keyGroup, valueGroup;

      HEMELB_MPI_CALL(MPI_Comm_group, (*commPtr, &keyGroup));
      HEMELB_MPI_CALL(MPI_Comm_group, (valueComm, &valueGroup));

      auto const N = Size();
      std::vector<int> keys(N), values(N);
      std::iota(keys.begin(), keys.end(), 0);

      HEMELB_MPI_CALL(MPI_Group_translate_ranks,
		      (keyGroup, N, keys.data(), valueGroup, values.data()));
      std::map<int, int> result;
      for (auto const & item : util::zip(keys, values))
      {
        result[std::get<0>(item)] = std::get<1>(item);
      }

      HEMELB_MPI_CALL(MPI_Group_free, (&keyGroup));
      HEMELB_MPI_CALL(MPI_Group_free, (&valueGroup));

      return result;
    }

}
