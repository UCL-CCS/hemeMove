// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_COLLISIONS_NORMAL_H
#define HEMELB_LB_COLLISIONS_NORMAL_H

#include "lb/collisions/BaseCollision.h"
#include "lb/kernels/BaseKernel.h"

namespace hemelb::lb
{
      /**
       * Normal collisions - use a formula to relax the distribution towards equilibrium.
       */
      template<typename KernelType>
      class Normal : public BaseCollision<Normal<KernelType>, KernelType>
      {
        public:
          using CKernel = KernelType;

          Normal(InitParams& initParams) :
              kernel(initParams)
          {
          }

          inline void DoCalculatePreCollision(HydroVars<KernelType>& hydroVars,
                                              const geometry::Site<geometry::Domain>& site)
          {
            kernel.CalculateDensityMomentumFeq(hydroVars, site.GetIndex());
          }

          inline void DoCollide(const LbmParameters* lbmParams,
                                HydroVars<KernelType>& iHydroVars)
          {
            kernel.Collide(lbmParams, iHydroVars);
          }

          KernelType kernel;

      };

}

#endif /* HEMELB_LB_COLLISIONS_NORMAL_H */
