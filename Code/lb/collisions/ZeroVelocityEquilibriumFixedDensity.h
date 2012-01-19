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
              BaseCollision<ZeroVelocityEquilibriumFixedDensity<KernelType>, KernelType>(), kernel(initParams), boundaryObject(initParams.boundaryObject)
          {
          }

          inline void DoCalculatePreCollision(kernels::HydroVars<KernelType>& hydroVars,
                                              const geometry::Site& site)
          {
            hydroVars.density = boundaryObject->GetBoundaryDensity(site.GetBoundaryId());

            hydroVars.v_x = 0.0;
            hydroVars.v_y = 0.0;
            hydroVars.v_z = 0.0;

            kernel.CalculateFeq(hydroVars, site.GetIndex());
          }

          inline void DoCollide(const LbmParameters* lbmParams,
                                kernels::HydroVars<KernelType>& iHydroVars)
          {
            for (Direction direction = 0; direction < D3Q15::NUMVECTORS; ++direction)
            {
              iHydroVars.GetFPostCollision()[direction] = iHydroVars.GetFEq()[direction];
            }
          }

          inline void DoReset(kernels::InitParams* init)
          {
            kernel.Reset(init);
          }

        private:
          KernelType kernel;
          boundaries::BoundaryValues* boundaryObject;

      };

    }
  }
}

#endif /* HEMELB_LB_COLLISIONS_ZEROVELOCITYEQUILIBRIUMFIXEDDENSITY_H */
