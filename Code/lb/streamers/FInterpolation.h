#ifndef HEMELB_LB_STREAMERS_FINTERPOLATION_H
#define HEMELB_LB_STREAMERS_FINTERPOLATION_H

#include "lb/streamers/BaseStreamer.h"

/**
 * TODO This class hasn't been modified to fit the new kernel / collision / streamer
 * hierarchy. This needs to happen before it can be used in the LBM.
 */

namespace hemelb
{
  namespace lb
  {
    namespace streamers
    {

      template<typename tCollisionOperator>
      class FInterpolation : public Implementation
      {

        public:
          template<bool tDoRayTracing>
          void DoStreamAndCollide(WallCollision* mWallCollision,
                                  const site_t iFirstIndex,
                                  const site_t iSiteCount,
                                  const LbmParameters* iLbmParams,
                                  geometry::LatticeData* bLatDat,
                                  hemelb::vis::Control *iControl);

          template<bool tDoRayTracing>
          void DoPostStep(WallCollision* mWallCollision,
                          const site_t iFirstIndex,
                          const site_t iSiteCount,
                          const LbmParameters* iLbmParams,
                          geometry::LatticeData* bLatDat,
                          hemelb::vis::Control *iControl);

      };

      template<typename tCollisionOperator>
      template<bool tDoRayTracing>
      void FInterpolation<tCollisionOperator>::DoStreamAndCollide(WallCollision* mWallCollision,
                                                                  const site_t iFirstIndex,
                                                                  const site_t iSiteCount,
                                                                  const LbmParameters* iLbmParams,
                                                                  geometry::LatticeData* bLatDat,
                                                                  hemelb::vis::Control *iControl)
      {
        for (site_t lIndex = iFirstIndex; lIndex < (iFirstIndex + iSiteCount); lIndex++)
        {
          distribn_t* f = bLatDat->GetFOld(lIndex * D3Q15::NUMVECTORS);
          distribn_t density, v_x, v_y, v_z, f_neq[15];
          // Temporarily store f_eq in f_neq. Rectified later.
          D3Q15::CalculateDensityVelocityFEq(f, density, v_x, v_y, v_z, f_neq);

          for (unsigned int ii = 0; ii < D3Q15::NUMVECTORS; ii++)
          {
            * (bLatDat->GetFNew(bLatDat->GetStreamedIndex(lIndex, ii))) = f[ii] += iLbmParams->Omega
                * (f_neq[ii] = f[ii] - f_neq[ii]);
          }

          UpdateMinsAndMaxes < tDoRayTracing
              > (v_x, v_y, v_z, lIndex, f_neq, density, bLatDat, iLbmParams, iControl);
        }
      }

      template<typename tCollisionOperator>
      template<bool tDoRayTracing>
      void FInterpolation<tCollisionOperator>::DoPostStep(WallCollision* mWallCollision,
                                                          const site_t iFirstIndex,
                                                          const site_t iSiteCount,
                                                          const LbmParameters* iLbmParams,
                                                          geometry::LatticeData* bLatDat,
                                                          hemelb::vis::Control *iControl)
      {
        for (site_t lIndex = iFirstIndex; lIndex < (iFirstIndex + iSiteCount); lIndex++)
        {
          // Handle odd indices, then evens - it's slightly easier to take the odd
          // and even cases separately.
          for (unsigned int l = 1; l < D3Q15::NUMVECTORS; l++)
          {
            if (bLatDat->HasBoundary(lIndex, l))
            {
              double twoQ = 2.0 * bLatDat->GetCutDistance(lIndex, l);

              * (bLatDat->GetFNew(lIndex * D3Q15::NUMVECTORS + D3Q15::INVERSEDIRECTIONS[l])) =
                  (twoQ < 1.0)
                  ?
                    (*bLatDat->GetFNew(lIndex * D3Q15::NUMVECTORS + l)
                        + twoQ
                            * (*bLatDat->GetFOld(lIndex * D3Q15::NUMVECTORS + l)
                                - *bLatDat->GetFNew(lIndex * D3Q15::NUMVECTORS + l)))
                        :
                    (*bLatDat->GetFOld(lIndex * D3Q15::NUMVECTORS + D3Q15::INVERSEDIRECTIONS[l])
                        + (1. / twoQ)
                            * (*bLatDat->GetFOld(lIndex * D3Q15::NUMVECTORS + l)
                                - *bLatDat->GetFOld(lIndex * D3Q15::NUMVECTORS
                                    + D3Q15::INVERSEDIRECTIONS[l])));
            }
          }
        }
      }

    }
  }
}

#endif /* HEMELB_LB_STREAMERS_FINTERPOLATION_H */
