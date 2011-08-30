#ifndef HEMELB_LB_KERNELS_LBGK_H
#define HEMELB_LB_KERNELS_LBGK_H

#include <cstdlib>
#include "lb/HFunction.h"
#include "util/utilityFunctions.h"
#include "lb/kernels/BaseKernel.h"

namespace hemelb
{
  namespace lb
  {
    namespace kernels
    {
      /**
       * LBGK: This class implements the LBGK single-relaxation time kernel.
       */
      class LBGK : public BaseKernel<LBGK>
      {
        public:
          LBGK(InitParams& initParams)
          {
          }

          void DoCalculateDensityVelocityFeq(HydroVars<LBGK>& hydroVars, site_t index)
          {
            D3Q15::CalculateDensityVelocityFEq(hydroVars.f,
                                               hydroVars.density,
                                               hydroVars.v_x,
                                               hydroVars.v_y,
                                               hydroVars.v_z,
                                               hydroVars.f_eq);
          }

          void DoCalculateFeq(HydroVars<LBGK>& hydroVars, site_t index)
          {
            D3Q15::CalculateFeq(hydroVars.density,
                                hydroVars.v_x,
                                hydroVars.v_y,
                                hydroVars.v_z,
                                hydroVars.f_eq);
          }

          distribn_t DoCollide(const LbmParameters* const lbmParams,
                               HydroVars<LBGK>& hydroVars,
                               unsigned int direction)
          {
            return hydroVars.f[direction] + hydroVars.f_neq[direction] * lbmParams->Omega;
          }
      };

    }
  }
}

#endif /* HEMELB_LB_COLLISIONS_IMPLEMENTATIONS_LBGK_H */
