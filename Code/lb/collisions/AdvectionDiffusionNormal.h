
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_COLLISIONS_ADVECTIONDIFFUSIONNORMAL_H
#define HEMELB_LB_COLLISIONS_ADVECTIONDIFFUSIONNORMAL_H

#include "lb/collisions/BaseCollision.h"
#include "lb/kernels/BaseKernel.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {
      /**
       * Normal collisions - use a formula to relax the distribution towards equilibrium.
       */
      template<typename KernelType>
      class AdvectionDiffusionNormal : public AdvectionDiffusionBaseCollision<Normal<KernelType>, KernelType>
      {
        public:
          typedef KernelType CKernel;

          AdvectionDiffusionNormal(kernels::InitParams& initParams) :
              kernel(initParams)
          {
          }

          inline void DoCalculatePreCollision(kernels::HydroVars<KernelType>& hydroVars,
                                              lb::MacroscopicPropertyCache& propertyCache,
                                              const geometry::Site<geometry::LatticeData>& site)
          {
            kernel.CalculateDensityMomentumFeq(hydroVars, propertyCache, site.GetIndex());
          }

          inline void DoCollide(const LbmParameters* lbmParams,
                                kernels::HydroVars<KernelType>& iHydroVars)
          {
            kernel.Collide(lbmParams, iHydroVars);
          }


          KernelType kernel;

      };

    }
  }
}

#endif /* HEMELB_LB_COLLISIONS_ADVECTIONDIFFUSIONNORMAL_H */
