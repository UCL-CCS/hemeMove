#ifndef HEMELB_LB_COLLISIONS_IMPLZEROVELOCITYEQUILIBRIUM_H
#define HEMELB_LB_COLLISIONS_IMPLZEROVELOCITYEQUILIBRIUM_H

#include "lb/collisions/WallCollision.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {

      class ImplZeroVelocityEquilibrium : public WallCollision
      {
        public:
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
      };

    }
  }
}
#endif /* HEMELB_LB_COLLISIONS_IMPLZEROVELOCITYEQUILIBRIUM_H */
