// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_STREAMERS_STREAMERTYPEFACTORY_H
#define HEMELB_LB_STREAMERS_STREAMERTYPEFACTORY_H

#include "lb/streamers/Common.h"
#include "lb/streamers/BulkStreamer.h"

namespace hemelb::lb
{
    /**
     * Template to produce Streamers that can cope with fluid-fluid,
     * fluid-wall and fluid-iolet links. Requires two classes as arguments:
     * 1) a link streamer that will handle the wall links, and
     * 2) a link streamer that will handle the iolet links.
     * In the case that one of these isn't needed, supply NullLink<Collision>.
     * Note that both link streamers must have the same CollisionType.
     */
    template<link_streamer WallLinkImpl, link_streamer IoletLinkImpl>
    requires std::same_as<typename WallLinkImpl::CollisionType, typename IoletLinkImpl::CollisionType>
    class StreamerTypeFactory
    {
    public:
        using CollisionType = typename WallLinkImpl::CollisionType;
        using VarsType = typename CollisionType::VarsType;
        using LatticeType = typename CollisionType::LatticeType;

    private:
        // Use these in the if statements below so the compiler can optimise them away if false.
        static constexpr bool can_have_wall = !std::same_as<WallLinkImpl, NullLink<CollisionType>>;
        static constexpr bool can_have_iolet = !std::same_as<IoletLinkImpl, NullLink<CollisionType>>;

        CollisionType collider;
        BulkLink<CollisionType> bulkLinkDelegate;
        WallLinkImpl wallLinkDelegate;
        IoletLinkImpl ioletLinkDelegate;

    public:
        StreamerTypeFactory(InitParams& initParams) :
                collider(initParams), bulkLinkDelegate(collider, initParams),
                wallLinkDelegate(collider, initParams), ioletLinkDelegate(collider, initParams)
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

                for (Direction ii = 0; ii < LatticeType::NUMVECTORS; ii++)
                {
                    if (can_have_iolet && site.HasIolet(ii))
                    {
                        ioletLinkDelegate.StreamLink(lbmParams, latDat, site, hydroVars, ii);
                    }
                    else if (can_have_wall && site.HasWall(ii))
                    {
                        wallLinkDelegate.StreamLink(lbmParams, latDat, site, hydroVars, ii);
                    }
                    else
                    {
                        bulkLinkDelegate.StreamLink(lbmParams, latDat, site, hydroVars, ii);
                    }
                }

                UpdateCachePostCollision(site,
                                         hydroVars,
                                         lbmParams,
                                         propertyCache);
            }
        }

        void PostStep(const site_t beginIndex, const site_t endIndex,
                      const LbmParameters* lbmParams, geometry::FieldData& latticeData,
                      lb::MacroscopicPropertyCache& propertyCache)
        {
            for (site_t siteIndex = beginIndex; siteIndex < endIndex; ++siteIndex)
            {
                geometry::Site<geometry::FieldData> site = latticeData.GetSite(siteIndex);
                for (unsigned int direction = 0; direction < LatticeType::NUMVECTORS; direction++)
                {
                    if (can_have_wall && site.HasWall(direction))
                    {
                        wallLinkDelegate.PostStepLink(latticeData, site, direction);
                    }
                    else if (can_have_iolet && site.HasIolet(direction))
                    {
                        ioletLinkDelegate.PostStepLink(latticeData, site, direction);
                    }
                }
            }

        }
    };
}
#endif
