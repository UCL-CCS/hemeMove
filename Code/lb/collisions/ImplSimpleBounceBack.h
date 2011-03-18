#ifndef HEMELB_LB_COLLISIONS_IMPLSIMPLEBOUNCEBACK_H
#define HEMELB_LB_COLLISIONS_IMPLSIMPLEBOUNCEBACK_H

#include "lb/collisions/WallCollision.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {

      class ImplSimpleBounceBack : public WallCollision
      {
        public:
          void DoCollisions(const bool iDoRayTracing,
                            const int iFirstIndex,
                            const int iSiteCount,
                            const LbmParameters &iLbmParams,
                            MinsAndMaxes &bMinimaAndMaxima,
                            geometry::LocalLatticeData &bLocalLatDat,
                            hemelb::vis::Control *iControl);

        private:
          template<bool tDoRayTracing>
          void DoCollisionsInternal(const int iFirstIndex,
                                    const int iSiteCount,
                                    const LbmParameters &iLbmParams,
                                    MinsAndMaxes &bMinimaAndMaxima,
                                    geometry::LocalLatticeData &bLocalLatDat,
                                    hemelb::vis::Control *iControl);

      };

    }
  }
}
#endif /* HEMELB_LB_COLLISIONS_IMPLSIMPLEBOUNCEBACK_H */
