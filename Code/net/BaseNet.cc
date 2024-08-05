// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

/*! \file net.cc
 \brief In this file the functions useful to discover the topology used and
 to create and delete the domain decomposition and the various
 buffers are defined.
 */

#include <cstdlib>
#include <cmath>
#include <cstdio>

#include "net/BaseNet.h"
#include "util/numerical.h"
#include "util/Vector3D.h"
#include "net/IOCommunicator.h"
namespace hemelb
{
  namespace net
  {

    void BaseNet::Dispatch()
    {
      Send();
      Receive();
      Wait();
    }

    BaseNet::BaseNet(const MpiCommunicator &commObject) :
        BytesSent(0), SyncPointsCounted(0), communicator(commObject)
    {
    }

    void BaseNet::Receive()
    {
      ReceiveGathers();
      ReceiveGatherVs();
      ReceiveAllToAll();
      // Ensure collectives are called before point-to-point, as some implementing mixins implement collectives via point-to-point
      ReceivePointToPoint();
    }

    void BaseNet::Send()
    {

      SendGathers();
      SendGatherVs();
      SendAllToAll();
      // Ensure collectives are called before point-to-point, as some implementing mixins implement collectives via point-to-point
      SendPointToPoint();
    }

    void BaseNet::Wait()
    {
      SyncPointsCounted++; //DTMP: counter for monitoring purposes.

      WaitGathers();
      WaitGatherVs();
      WaitPointToPoint();
      WaitAllToAll();

      displacementsBuffer.clear();
      countsBuffer.clear();
    }

    std::vector<int> & BaseNet::GetDisplacementsBuffer()
    {
      displacementsBuffer.push_back(std::vector<int>());
      return displacementsBuffer.back();
    }

    std::vector<int> & BaseNet::GetCountsBuffer()
    {
      countsBuffer.push_back(std::vector<int>());
      return countsBuffer.back();
    }

  }
}
