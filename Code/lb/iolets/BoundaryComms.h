// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_IOLETS_BOUNDARYCOMMS_H
#define HEMELB_LB_IOLETS_BOUNDARYCOMMS_H

#include "geometry/Domain.h"
#include "lb/SimulationState.h"
#include "util/UnitConverter.h"
#include "lb/iolets/BoundaryCommunicator.h"

namespace hemelb::lb
{

      class BoundaryComms
      {
      public:
          BoundaryComms(SimulationState* iSimState, std::vector<int> &iProcsList,
                        const BoundaryCommunicator& boundaryComm, bool iHasBoundary);

          void Wait();

          // It is up to the caller to make sure only BCproc calls send
          void Send(distribn_t* density);
          void Receive(distribn_t* density);

          const std::vector<int>& GetListOfProcs() const
          {
            return procsList;
          }

          void ReceiveDoubles(double* double_array, int size);
          void WaitAllComms();
          void FinishSend();

      private:
          // MPI tag for communication
          static constexpr int BC_TAG = 100;
          // This is necessary to support BC proc having fluid sites
          bool hasBoundary;

          // These are only assigned on the BC proc as it is the only one that needs to know
          // which proc has which IOlet
          int nProcs;
          std::vector<int> procsList;
          const BoundaryCommunicator& bcComm;

          std::vector<net::MpiRequest> sendRequest;
          net::MpiRequest receiveRequest;

          SimulationState* mState;
      };

}

#endif /* HEMELB_LB_IOLETS_BOUNDARYCOMMS_H */
