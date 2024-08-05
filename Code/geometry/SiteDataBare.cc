// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "constants.h"
#include "geometry/SiteDataBare.h"
#include "Exception.h"

namespace hemelb::geometry
{
    /**
     * Constructor for a SiteData object
     *
     * NOTE: that this assumes a single fluid site will be next to at most one inlet or outlet.
     * This is a safe assumption in simple geometries but may not always be the case.
     *
     * @param readResult
     */
    SiteData::SiteData(const GeometrySite& readResult)
    {
      if (!readResult.isFluid)
      {
        type = SOLID_TYPE;
        wallIntersection = 0;
        ioletIntersection = 0;
        ioletId = -1;
      }
      else
      {
        ioletId = -1;
        wallIntersection = 0;
        ioletIntersection = 0;

        bool hadInlet = false;
        bool hadOutlet = false;

        // Iterate over each (non-zero) direction
        for (Direction direction = 1; direction <= readResult.links.size(); ++direction)
        {
	  using CutType = io::formats::geometry::CutType;
          // Get the link
          const GeometrySiteLink& link = readResult.links[direction - 1];
	  switch (link.type) {
	  case CutType::WALL:
            wallIntersection |= 1 << (direction - 1);
	    break;

          case CutType::INLET:
	    // This is meant to fall-through.
	  case CutType::OUTLET:
            ioletId = link.ioletId;
            ioletIntersection |= 1 << (direction - 1);
	    if (link.type == CutType::INLET) {
	      hadInlet = true;
	    } else {
	      hadOutlet = true;
	    }
	    break;
	  default:
	    break;
	  }
        }

        type = hadInlet
          ? INLET_TYPE
          : (hadOutlet
            ? OUTLET_TYPE
            : FLUID_TYPE);

      }
    }

    bool SiteData::IsWall() const
    {
      return wallIntersection != 0;
    }

    bool SiteData::IsSolid() const
    {
      return GetSiteType() == SOLID_TYPE;
    }

    unsigned SiteData::GetCollisionType() const
    {
      if (wallIntersection == 0)
      {
        // No solid wall intersections
        switch (type)
        {
          case FLUID_TYPE:
            return FLUID;

          case INLET_TYPE:
            return INLET;

          case OUTLET_TYPE:
            return OUTLET;

          case SOLID_TYPE:
            throw (Exception() << "Requesting collision type for solid site!");
        }
      }
      else
      {
        // There are solid wall intersections
        switch (type)
        {
          case FLUID_TYPE:
            return WALL;

          case INLET_TYPE:
            return INLET | WALL;

          case OUTLET_TYPE:
            return OUTLET | WALL;

          case SOLID_TYPE:
            throw (Exception() << "Requesting collision type for solid site!");

        }
      }

      throw (Exception() << "Requesting collision type for solid site!");
      // The end of this function should never be reached. Adding return statement to please CRAY compiler
      return 0u;
    }

    bool SiteData::HasWall(Direction direction) const
    {
      // If at the zero direction then always false, so set mask to all zeros.
      const unsigned mask = direction ?
	1U << (direction - 1U) :
	0U;
      return (wallIntersection & mask) != 0;
    }

    bool SiteData::HasIolet(Direction direction) const
    {
      // If at the zero direction then always false, so set mask to all zeros.
      unsigned mask = direction ?
	1U << (direction - 1U) :
	0U;
      return (ioletIntersection & mask) != 0;
    }

    uint32_t SiteData::GetIoletIntersectionData() const
    {
      return ioletIntersection;
    }

    uint32_t SiteData::GetWallIntersectionData() const
    {
      return wallIntersection;
    }

}
