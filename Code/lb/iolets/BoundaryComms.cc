// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "lb/iolets/BoundaryComms.h"
#include "lb/iolets/BoundaryValues.h"
#include "net/IOCommunicator.h"
#include "util/numerical.h"

namespace hemelb::lb
{

      BoundaryComms::BoundaryComms(SimulationState* iSimState, std::vector<int> &iProcsList,
                                   const BoundaryCommunicator& boundaryComm, bool iHasBoundary) :
          hasBoundary(iHasBoundary), nProcs((int) iProcsList.size()), procsList(iProcsList),
              bcComm(boundaryComm)
      {
        /* iProcsList contains the procs containing said Boundary/iolet, but NOT proc 0 (the BoundaryControlling/BC proc)! */

        // Only BC proc sends
        if (bcComm.IsCurrentProcTheBCProc())
          sendRequest.resize(nProcs);
      }

      void BoundaryComms::Wait()
      {
        if (hasBoundary)
            receiveRequest.Wait();
      }

      void BoundaryComms::WaitAllComms()
      {
        // Now wait for all to complete
        if (bcComm.IsCurrentProcTheBCProc()) {
          net::MpiRequest::Waitall(sendRequest);

          if (hasBoundary)
              receiveRequest.Wait();
        } else {
            receiveRequest.Wait();
        }

      }

      // It is up to the caller to make sure only BCproc calls send
      void BoundaryComms::Send(distribn_t* density)
      {
        for (int proc = 0; proc < nProcs; proc++)
            sendRequest[proc] = bcComm.Issend(*density, procsList[proc], BC_TAG);
      }

      void BoundaryComms::Receive(distribn_t* density)
      {
        if (hasBoundary)
            receiveRequest = bcComm.Irecv(*density, bcComm.GetBCProcRank(), BC_TAG);
      }

      void BoundaryComms::FinishSend()
      {
        // Don't move on to next step with BC proc until all messages have been sent
        // Precautionary measure to make sure proc doesn't overwrite, before message is sent
        if (bcComm.IsCurrentProcTheBCProc())
            net::MpiRequest::Waitall(sendRequest);
      }

}
