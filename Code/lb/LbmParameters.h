#ifndef HEMELB_LB_LBMPARAMETERS_H
#define HEMELB_LB_LBMPARAMETERS_H

#include <cmath>
#include "constants.h"

namespace hemelb
{
  namespace lb
  {
    enum StressTypes
    {
      VonMises = 0,
      ShearStress = 1,
      IgnoreStress = 2
    };

    struct LbmParameters
    {
      public:
        LbmParameters(distribn_t timeStepLength, distribn_t voxelSize)
        {
          Update(timeStepLength, voxelSize);
        }

        void Update(distribn_t timeStepLength, distribn_t voxelSize)
        {
          timestep = timeStepLength;
          tau = 0.5 + (timeStepLength * BLOOD_VISCOSITY_Pa_s / BLOOD_DENSITY_Kg_per_m3) / (Cs2
              * voxelSize * voxelSize);

          omega = -1.0 / tau;
          stressParameter = (1.0 - 1.0 / (2.0 * tau)) / sqrt(2.0);
          beta = -1.0 / (2.0 * tau);
        }

        // TODO naming convention question: shall this method be called TimeStep like the rest in the struct or start with a verb?
        distribn_t GetTimeStep() const
        {
          return timestep;
        }

        distribn_t Omega() const
        {
          return omega;
        }

        distribn_t Tau() const
        {
          return tau;
        }

        distribn_t StressParameter() const
        {
          return stressParameter;
        }

        distribn_t Beta() const
        {
          return beta;
        }

        StressTypes StressType;

      private:
        distribn_t timestep;
        distribn_t omega;
        distribn_t tau;
        distribn_t stressParameter;
        distribn_t beta; // Viscous dissipation in ELBM
    };
  }
}

#endif //HEMELB_LB_LBMPARAMETERS_H
