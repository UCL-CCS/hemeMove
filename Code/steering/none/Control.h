#ifndef HEMELB_STEERING_NONE_CONTROL_H
#define HEMELB_STEERING_NONE_CONTROL_H

#include "lb/SimulationState.h"
#include "lb/lb.h"
#include "vis/Control.h"

namespace hemelb
{
  namespace steering
  {

    class Control
    {
      public:
        // Singleton
        static Control* Init(bool isCurrentProcTheSteeringProc);
        static Control* Get(void);

        void StartNetworkThread(lb::LBM* lbm,
                                lb::SimulationState *iSimState,
                                const lb::LbmParameters *iLbmParams);

        void
            UpdateSteerableParameters(bool shouldRenderForSnapshot,
                                      int* perform_rendering,
                                      hemelb::lb::SimulationState &iSimulationState,
                                      hemelb::vis::Control* visController,
                                      lb::LBM* lbm);
        bool ShouldRenderForNetwork();

      protected:
        // Singleton pattern
        static Control* singleton;

        Control(bool isCurrentProcTheSteeringProc);
        ~Control();
    };

  }

}

#endif /* HEMELB_STEERING_NONE_CONTROL_H */
