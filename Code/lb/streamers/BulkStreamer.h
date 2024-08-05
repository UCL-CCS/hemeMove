// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_STREAMERS_BULKSTREAMER_H
#define HEMELB_LB_STREAMERS_BULKSTREAMER_H

#include "geometry/FieldData.h"
#include "lb/streamers/Common.h"
#include "lb/HFunction.h"

namespace hemelb::lb
{
    /// Do bulk streaming along a link
    /// Instantiations satisfy link_streamer
    template<collision_type C>
    class BulkLink
    {
    public:
        using CollisionType = C;
        using KernelType = typename CollisionType::KernelType;
        using VarsType = typename CollisionType::VarsType;
        using LatticeType = typename KernelType::LatticeType;

        BulkLink(CollisionType& delegatorCollider,
                 InitParams& initParams)
        {
        }

        void StreamLink(const LbmParameters* lbmParams,
                        geometry::FieldData& latticeData,
                        const geometry::Site<geometry::FieldData>& site,
                        VarsType& hydroVars,
                        const Direction& direction)
        {
            * (latticeData.GetFNew(site.GetStreamedIndex<LatticeType>(direction))) =
                    hydroVars.GetFPostCollision()[direction];
        }

        void PostStepLink(geometry::FieldData& latticeData,
                          const geometry::Site<geometry::FieldData>& site,
                          const Direction& direction) {
        }
    };

    template<collision_type C>
    class BulkStreamer
    {
    public:
        using CollisionType = C;
        using KernelType = typename C::KernelType;
        using LatticeType = typename C::LatticeType;
        using VarsType = typename C::VarsType;

    private:
        CollisionType collider;
        BulkLink<CollisionType> bulkLinkDelegate;
        static_assert(link_streamer<BulkLink<CollisionType>>);

    public:
        BulkStreamer(InitParams& initParams) :
                collider(initParams), bulkLinkDelegate(collider, initParams)
        {
        }

        void StreamAndCollide(const site_t beginIndex, const site_t endIndex,
                              const LbmParameters* lbmParams,
                              geometry::FieldData& latDat,
                              lb::MacroscopicPropertyCache& propertyCache)
        {
            for (site_t siteIndex = beginIndex; siteIndex < endIndex; ++siteIndex)
            {
                geometry::Site<geometry::FieldData> site = latDat.GetSite(siteIndex);
                VarsType hydroVars(site);

                ///< @todo #126 This value of tau will be updated by some kernels within the collider code (e.g. LBGKNN). It would be nicer if tau is handled in a single place.
                hydroVars.tau = lbmParams->GetTau();

                collider.CalculatePreCollision(hydroVars, site);

                collider.Collide(lbmParams, hydroVars);

                for (unsigned int ii = 0; ii < LatticeType::NUMVECTORS; ii++)
                {
                    bulkLinkDelegate.StreamLink(lbmParams, latDat, site, hydroVars, ii);
                }

                UpdateCachePostCollision(site, hydroVars, lbmParams, propertyCache);
            }
        }

        void PostStep(const site_t beginIndex, const site_t iSiteCount,
                      const LbmParameters* iLbmParams, geometry::FieldData& bLatDat,
                      lb::MacroscopicPropertyCache& propertyCache)
        {
        }

    };
}
#endif
