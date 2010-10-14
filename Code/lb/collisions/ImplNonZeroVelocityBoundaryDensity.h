#ifndef HEMELB_LB_COLLISIONS_IMPLNONZEROVELOCITYBOUNDARYDENSITY_H
#define HEMELB_LB_COLLISIONS_IMPLNONZEROVELOCITYBOUNDARYDENSITY_H

#include "lb/collisions/InletOutletCollision.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {

      class ImplNonZeroVelocityBoundaryDensity : public InletOutletCollision
      {
        public:
          ImplNonZeroVelocityBoundaryDensity(double* iBounaryDensityArray);

          void DoCollisions(double omega, int i, double *density, double *v_x,
            double *v_y, double *v_z, double f_neq[], Net* net);

        private:
          double* mBoundaryDensityArray;
      };

    }
  }
}

#endif /* HEMELB_LB_COLLISIONS_IMPLNONZEROVELOCITYBOUNDARYDENSITY_H */
