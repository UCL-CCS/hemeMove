// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_GEOMETRY_SITEDATABARE_H
#define HEMELB_GEOMETRY_SITEDATABARE_H

#include "units.h"
#include "geometry/GeometrySite.h"
#include "geometry/SiteType.h"

namespace hemelb::geometry
{
    class SiteData
    {
    public:
        SiteData(const GeometrySite& siteReadResult);
        // Default constructor allows one to use operator[] for maps
        SiteData() = default;

        bool IsWall() const;
        bool IsSolid() const;
        unsigned GetCollisionType() const;

        inline SiteType GetSiteType() const
        {
          return type;
        }

        inline SiteType& GetSiteType()
        {
          return type;
        }

        inline int GetIoletId() const
        {
          return ioletId;
        }

        inline int& GetIoletId()
        {
          return ioletId;
        }

        bool HasWall(Direction direction) const;
        bool HasIolet(Direction direction) const;

        /**
         * These functions return internal representations and should only be used for debugging.
         */
        uint32_t GetIoletIntersectionData() const;
        inline uint32_t &GetIoletIntersectionData()
        {
          return ioletIntersection;
        }
        uint32_t GetWallIntersectionData() const;
        inline uint32_t &GetWallIntersectionData()
        {
          return wallIntersection;
        }

    protected:
        /**
         * This is a bit mask for whether a wall is hit by links in each direction.
         */
        uint32_t wallIntersection = 0U;

        /**
         * This is a bit mask for whether an iolet is hit by links in each direction.
         */
        uint32_t ioletIntersection = 0U;

        SiteType type = SOLID_TYPE;
        int ioletId = -1;
    };
}

#endif /* HEMELB_GEOMETRY_SITEDATABARE_H */
