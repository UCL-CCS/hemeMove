#include <limits>
#include "lb/collisions/Collision.h"

namespace hemelb
{
  namespace lb
  {
    namespace collisions
    {

      // Default constructor, implemented so we can make it protected
      // and prevent instantiation of this base class.
      Collision::Collision()
      {

      }
      Collision::~Collision()
      {
      }

      void Collision::DoCollisions(const bool iDoRayTracing,
                                   const site_t iFirstIndex,
                                   const site_t iSiteCount,
                                   const LbmParameters* iLbmParams,
                                   geometry::LatticeData* bLatDat,
                                   hemelb::vis::Control *iControl)
      {
        // Standard implementation - do nothing.
      }

      void Collision::PostStep(const bool iDoRayTracing,
                               const site_t iFirstIndex,
                               const site_t iSiteCount,
                               const LbmParameters* iLbmParams,
                               geometry::LatticeData* bLatDat,
                               hemelb::vis::Control *iControl)
      {
        // Standard implementation - do nothing.
      }
    }
  }
}
