#ifndef HEMELB_GEOMETRY_SITE_H
#define HEMELB_GEOMETRY_SITE_H

#include "units.h"
#include "geometry/SiteData.h"
#include "util/Vector3D.h"
#include "util/ConstSelector.h"

namespace hemelb
{
  namespace geometry
  {
    // Forward definition of LatticeData; makes things easier.
    class LatticeData;

    template<bool isConst>
    class BaseSite
    {
      public:
        BaseSite(site_t localContiguousIndex, typename util::constSelector<isConst, geometry::LatticeData&,
            const geometry::LatticeData&>::type latticeData) :
          index(localContiguousIndex), latticeData(latticeData)
        {
        }

        inline bool IsEdge() const
        {
          return GetSiteData().IsEdge();
        }

        inline bool IsSolid() const
        {
          return GetSiteData().IsSolid();
        }

        inline unsigned GetCollisionType() const
        {
          return GetSiteData().GetCollisionType();
        }

        inline SiteType GetSiteType() const
        {
          return GetSiteData().GetSiteType();
        }

        inline int GetBoundaryId() const
        {
          return GetSiteData().GetBoundaryId();
        }

        inline bool HasBoundary(Direction direction) const
        {
          return GetSiteData().HasBoundary(direction);
        }

        template<typename LatticeType>
        inline distribn_t GetWallDistance(Direction direction) const
        {
          return latticeData.template GetCutDistance<LatticeType> (index, direction);
        }

        inline const util::Vector3D<distribn_t>& GetWallNormal() const
        {
          return latticeData.GetNormalToWall(index);
        }

        inline site_t GetIndex() const
        {
          return index;
        }

        /**
         * This returns the index of the distribution to stream to. If streaming would take the
         * distribution out of the geometry, we instead stream to the 'rubbish site', an extra
         * position in the array that doesn't correspond to any site in the geometry.
         *
         * @param direction
         * @return
         */
        template<typename LatticeType>
        inline site_t GetStreamedIndex(Direction direction) const
        {
          return latticeData.template GetStreamedIndex<LatticeType> (index, direction);
        }

        template<typename LatticeType>
        inline typename util::constSelector<isConst, distribn_t*, const distribn_t*>::type GetFOld() const
        {
          return latticeData.GetFOld(index * LatticeType::NUMVECTORS);
        }

        inline const SiteData GetSiteData() const
        {
          return latticeData.GetSiteData(index);
        }

        inline const util::Vector3D<site_t>& GetGlobalSiteCoords() const
        {
          return latticeData.GetGlobalSiteCoords(index);
        }

      private:
        site_t index;
        typename util::constSelector<isConst, LatticeData&, const LatticeData&>::type latticeData;
    };

    typedef BaseSite<false> Site;

    typedef BaseSite<true> ConstSite;
  }
}

#endif
