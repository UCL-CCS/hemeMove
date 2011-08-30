#ifndef HEMELB_LB_COLLISIONS_ZEROVELOCITYEQUILIBRIUMFIXEDDENSITY_H
#define HEMELB_LB_COLLISIONS_ZEROVELOCITYEQUILIBRIUMFIXEDDENSITY_H

#include "lb/collisions/BaseCollision.h"
#include "lb/kernels/BaseKernel.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {
      /**
       * Collision operator that uses an externally imposed density, zero velocity and relaxes
       * fully to equilibrium.
       */
      template<typename KernelType>
      class ZeroVelocityEquilibriumFixedDensity : public BaseCollision<
          ZeroVelocityEquilibriumFixedDensity<KernelType>, KernelType>
      {
        public:
          typedef KernelType CKernel;

          ZeroVelocityEquilibriumFixedDensity(kernels::InitParams& initParams) :
              BaseCollision<ZeroVelocityEquilibriumFixedDensity<KernelType>, KernelType>(), kernel(initParams), boundaryDensityArray(initParams.boundaryDensityArray), latDat(initParams.latDat)
          {
          }

          void DoCalculatePreCollision(kernels::HydroVars<KernelType>& hydroVars,
                                       const site_t index)
          {
            hydroVars.density = boundaryDensityArray[latDat->GetBoundaryId(index)];

            hydroVars.v_x = 0.0;
            hydroVars.v_y = 0.0;
            hydroVars.v_z = 0.0;

            kernel.CalculateFeq(hydroVars, index);
          }

          distribn_t DoCollide(const LbmParameters* lbmParams,
                               unsigned int directionIndex,
                               kernels::HydroVars<KernelType>& iHydroVars)
          {
            return iHydroVars.f_eq[directionIndex];
          }

        private:
          KernelType kernel;
          distribn_t* boundaryDensityArray;
          const geometry::LatticeData* latDat;

      };

    }
  }
}

#endif /* HEMELB_LB_COLLISIONS_ZEROVELOCITYEQUILIBRIUMFIXEDDENSITY_H */
