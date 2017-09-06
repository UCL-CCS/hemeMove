
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_LBMPARAMETERS_H
#define HEMELB_LB_LBMPARAMETERS_H

#include <cmath>
#include "constants.h"
#include <vector>
#include <cassert>

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
        LbmParameters(PhysicalTime timeStepLength, PhysicalDistance voxelSize)
        {
          Update(timeStepLength, voxelSize);
        }

        void Update(PhysicalTime timeStepLength, PhysicalDistance voxelSizeMetres)
        {
          timestep = timeStepLength;
          voxelSize = voxelSizeMetres;
          tau = 0.5
              + (timeStepLength * BLOOD_VISCOSITY_Pa_s / BLOOD_DENSITY_Kg_per_m3)
                  / (Cs2 * voxelSize * voxelSize);

          omega = -1.0 / tau;
          advectionDiffusionTau = 0.5
              + (timeStepLength * DIFFUSIVITY_m2_per_s)
                  / (Cs2 * voxelSize * voxelSize);

          advectionDiffusionOmega = -1.0 / advectionDiffusionTau;
          stressParameter = (1.0 - 1.0 / (2.0 * tau)) / sqrt(2.0);
          beta = -1.0 / (2.0 * tau);
          alpha = KINETIC_RATE_m_per_s * (timeStepLength / voxelSize);
        }

        PhysicalTime GetTimeStep() const
        {
          return timestep;
        }

        PhysicalDistance GetVoxelSize() const
        {
          return voxelSize;
        }

        distribn_t GetOmega() const
        {
          return omega;
        }

        distribn_t GetTau() const
        {
          return tau;
        }

        distribn_t GetAdvectionDiffusionOmega() const
        {
          return advectionDiffusionOmega;
        }

        distribn_t GetAdvectionDiffusionTau() const
        {
          return advectionDiffusionTau;
        }

        distribn_t GetAdvectionDiffusionAlpha() const
        {
          return alpha;
        }

        distribn_t GetStressParameter() const
        {
          return stressParameter;
        }

        distribn_t GetBeta() const
        {
          return beta;
        }

        StressTypes StressType;

      private:
        PhysicalTime timestep;
        PhysicalDistance voxelSize;
        distribn_t omega;
        distribn_t tau;
        distribn_t alpha;
        distribn_t advectionDiffusionOmega;
        distribn_t advectionDiffusionTau;
        distribn_t stressParameter;
        distribn_t beta; ///< Viscous dissipation in ELBM
    };
  }
}

#endif //HEMELB_LB_LBMPARAMETERS_H
