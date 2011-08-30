#include "lb/streamers/InletOutletWallCollision.h"
#include "lb/collisions/CollisionVisitor.h"

namespace hemelb
{
  namespace lb
  {
    namespace streamers
    {
      InletOutletWallCollision::InletOutletWallCollision(boundaries::BoundaryValues* iBoundaryValues) :
        mBoundaryValues(iBoundaryValues)
      {

      }

      InletOutletWallCollision::~InletOutletWallCollision()
      {

      }

      distribn_t InletOutletWallCollision::getBoundaryDensityArray(const int index)
      {
        return mBoundaryValues->GetBoundaryDensity(index);
      }

      void InletOutletWallCollision::AcceptCollisionVisitor(collisions::CollisionVisitor* v,
                                                            const bool iDoRayTracing,
                                                            const site_t iFirstIndex,
                                                            const site_t iSiteCount,
                                                            const LbmParameters* iLbmParams,
                                                            geometry::LatticeData* bLatDat,
                                                            hemelb::vis::Control *iControl)
      {
        v->VisitInletOutletWall(this,
                                iDoRayTracing,
                                iFirstIndex,
                                iSiteCount,
                                iLbmParams,
                                bLatDat,
                                iControl);
      }

    }
  }
}
