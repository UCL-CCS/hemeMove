#ifndef HEMELB_STEERING_STEERINGCOMPONENT_H
#define HEMELB_STEERING_STEERINGCOMPONENT_H

#include "net/PhasedBroadcastRegular.h"
#include "lb/lb.h"
#include "lb/SimulationState.h"
#include "SimConfig.h"
#include "steering/Network.h"
#include "vis/DomainStats.h"
#include "vis/Control.h"

namespace hemelb
{
  namespace steering
  {
    namespace parameter
    {
      enum parameter 
      {
	SceneCentreX = 0,
	SceneCentreY = 1,
	SceneCentreZ = 2,
	Longitude = 3,
	Latitude = 4,
	Zoom = 5,
	Brightness = 6,
	PhysicalVelocityThresholdMax = 7,
	PhysicalStressThrehsholdMaximum = 8,
	PhysicalPressureThresholdMinimum = 9,
	PhysicalPressureThresholdMaximum = 10,
	GlyphLength = 11,
	PixelsX = 12,
	PixelsY = 13,
	NewMouseX = 14,
	NewMouseY = 15,
	SetIsTerminal = 16,
	Mode = 17,
	StreaklinePerPulsatilePeriod = 18,
	StreallineLength = 19,
	SetDoRendering = 20
      };
    }

    /**
     * SteeringComponent - class for passing steering data to all nodes.
     *
     * We pass this data at regular intervals. No initial action is required by all nodes, and
     * we only need to pass from the top-most node (which handles network communication) downwards,
     * on one iteration between each pair of consecutive depths.
     */
    class SteeringComponent : public net::PhasedBroadcastRegular<false, 1, 0, true, false>
    {
      public:
        SteeringComponent(Network* iNetwork,
                          vis::Control* iVisControl,
                          lb::LBM* iLbm,
                          net::Net * iNet,
                          lb::SimulationState * iSimState,
	                  SimConfig& iSimConfig);

        static bool RequiresSeparateSteeringCore();

        /*
         * This function initialises all of the steering parameters, on all nodes.
         */
        void Reset(SimConfig& iSimConfig);

        bool readyForNextImage;
        bool updatedMouseCoords;

      protected:
        void ProgressFromParent(unsigned long splayNumber);
        void ProgressToChildren(unsigned long splayNumber);

        void TopNodeAction();
        void Effect();

      private:
        void AssignValues();

        const static int STEERABLE_PARAMETERS = 20;
        const static unsigned int SPREADFACTOR = 10;

        bool isConnected;

        Network* mNetwork;
        lb::LBM* mLbm;
        lb::SimulationState* mSimState;
        vis::Control* mVisControl;
        float privateSteeringParams[STEERABLE_PARAMETERS + 1];
    };
  }
}

#endif /* HEMELB_STEERING_STEERINGCOMPONENT_H */
