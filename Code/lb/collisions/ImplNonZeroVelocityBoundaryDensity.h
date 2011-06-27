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
          ImplNonZeroVelocityBoundaryDensity(distribn_t* iBounaryDensityArray);

          void DoCollisions(const bool iDoRayTracing,
                            const site_t iFirstIndex,
                            const site_t iSiteCount,
                            const LbmParameters* iLbmParams,
                            geometry::LatticeData* bLatDat,
                            hemelb::vis::Control *iControl);

        private:
          template<bool tDoRayTracing>
          void DoCollisionsInternal(const site_t iFirstIndex,
                                    const site_t iSiteCount,
                                    const LbmParameters* iLbmParams,
                                    geometry::LatticeData* bLatDat,
                                    hemelb::vis::Control *iControl);

          distribn_t* mBoundaryDensityArray;

      };

    }
  }
}

#endif /* HEMELB_LB_COLLISIONS_IMPLNONZEROVELOCITYBOUNDARYDENSITY_H */
