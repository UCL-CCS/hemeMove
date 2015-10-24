//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_NET_MAPALLTOALL_H
#define HEMELB_NET_MAPALLTOALL_H

#include <map>
#include "net/mpi.h"

namespace hemelb
{
  namespace net
  {
    template <typename T>
    void MapAllToAll(const MpiCommunicator& comm,
                     const std::map<int, T>& valsToSend,
                     std::map<int, T>& receivedVals,
                     const int tag = 1234)
    {
      std::vector<MpiRequest> sendReqs(valsToSend.size());
      receivedVals.clear();
      int localRank = comm.Rank();

      typename std::map<int, T>::const_iterator iter = valsToSend.cbegin(),
          end = valsToSend.cend();
      for (int i = 0; iter != end; ++iter, ++i) {
        int rank = iter->first;
        const T& val = iter->second;
        if (rank == localRank)
        {
          receivedVals[localRank] = val;
        }
        else
        {
          // Synchronous to ensure we know when this is matched
          sendReqs[i] = comm.Issend(val, rank, tag);
        }
      }

      bool allProcsHaveHadAllSendsMatched = false;
      bool mySendsMatched = false;
      MpiRequest barrierReq;

      // Allocate outside of loop to ensure compiler isn't reinitialising this
      // all the time.
      MPI_Status status;
      while(!mySendsMatched || !allProcsHaveHadAllSendsMatched)
      {
        if (comm.Iprobe(MPI_ANY_SOURCE, tag, &status)) {
          // There is a message waiting
          comm.Recv(receivedVals[status.MPI_SOURCE], status.MPI_SOURCE, status.MPI_TAG);
        }

        if (!mySendsMatched)
        {
          if (MpiRequest::TestAll(sendReqs))
          {
            mySendsMatched = true;
            barrierReq = comm.Ibarrier();
          }
        }
        else
        {
          allProcsHaveHadAllSendsMatched = barrierReq.Test();
        }
      }
    }
  }
}

#endif /* HEMELB_NET_MAPALLTOALL_H */
