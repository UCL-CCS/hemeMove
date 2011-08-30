#ifndef HEMELB_LB_COLLISIONS_NONZEROVELOCITYEQUILIBRIUMFIXEDDENSITY_H
#define HEMELB_LB_COLLISIONS_NONZEROVELOCITYEQUILIBRIUMFIXEDDENSITY_H

#include "lb/collisions/BaseCollision.h"
#include "lb/kernels/BaseKernel.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {
      template<typename KernelType>
      class NonZeroVelocityEquilibriumFixedDensity : public BaseCollision<
          NonZeroVelocityEquilibriumFixedDensity<KernelType>, KernelType>
      {
        public:
          typedef KernelType CKernel;

          NonZeroVelocityEquilibriumFixedDensity(kernels::InitParams& initParams) :
              kernel(initParams), boundaryDensityArray(initParams.boundaryDensityArray), latDat(initParams.latDat)
          {

          }

          void DoCalculatePreCollision(kernels::HydroVars<KernelType>& hydroVars,
                                       const site_t index)
          {
            D3Q15::CalculateDensityAndVelocity(hydroVars.f,
                                               hydroVars.density,
                                               hydroVars.v_x,
                                               hydroVars.v_y,
                                               hydroVars.v_z);

            // Externally impose a density.
            hydroVars.density = boundaryDensityArray[latDat->GetBoundaryId(index)];

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

#endif /* HEMELB_LB_COLLISIONS_NONZEROVELOCITYEQUILIBRIUMFIXEDDENSITY_H */
