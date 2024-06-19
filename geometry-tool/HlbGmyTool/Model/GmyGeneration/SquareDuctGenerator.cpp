// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "SquareDuctGenerator.h"

#include "InconsistentFluidnessError.h"
#include "Neighbours.h"
#include "Site.h"

#include "Debug.h"

#include "io/formats/geometry.h"

#include <cmath>

namespace hemelb::gmytool::gmy {

using namespace io::formats;

SquareDuctGenerator::SquareDuctGenerator() : GeometryGenerator() {
  this->SquareDuct = new SquareDuctData;
}

SquareDuctGenerator::~SquareDuctGenerator() {
  delete this->SquareDuct;
}
void SquareDuctGenerator::ComputeBounds(double bounds[6]) const {
  for (int i = 0; i < 3; ++i) {
    bounds[2 * i] = this->SquareDuct->LowerBound[i];
    bounds[2 * i + 1] = this->SquareDuct->UpperBound[i];
  }
}

bool IsInsideDuct(SquareDuctData* duct, Vector& pt) {
  for (int i = 0; i < 3; ++i) {
    if ((pt[i] > duct->LowerBound[i]) && (pt[i] < duct->UpperBound[i])) {
      // OK
    } else {
      return false;
    }
  }
  return true;
}

/*
 * Given a site with known fluidness, examine the links to not-yet-visited
 * neighbouring sites. If the neighbours have unknown fluidness, set that.
 *
 * Since we wish to examine each link only once, this will set link properties
 * from neighbour => site as well as site => neighbour
 *
 */
void SquareDuctGenerator::ClassifySite(Site& site) {
  for (LaterNeighbourIterator neighIt = site.begin(); neighIt != site.end();
       ++neighIt) {
    Site& neigh = *neighIt;
    unsigned int iNeigh = neighIt.GetNeighbourIndex();
    int nHits;

    if (!neigh.IsFluidKnown) {
      neigh.IsFluidKnown = true;
      neigh.IsFluid = IsInsideDuct(this->SquareDuct, neigh.Position);
      if (neigh.IsFluid)
        neigh.CreateLinksVector();
    }

    // Four cases: fluid-fluid, solid-solid, fluid-solid and solid-fluid.
    // Will handle the last two together.
    if (site.IsFluid == neigh.IsFluid) {
      if (site.IsFluid) {
        // Fluid-fluid, must set CutType::NONE for both
        site.Links[iNeigh].Type = geometry::CutType::NONE;
        neigh.Links[Neighbours::inverses[iNeigh]].Type =
            geometry::CutType::NONE;
      } else {
        // solid-solid, nothing to do.
      }
    } else {
      // They differ, figure out which is fluid and which is solid.
      Site* fluid;
      Site* solid;

      // Index of the solid site from the fluid site.
      int iSolid;

      if (site.IsFluid) {
        fluid = &site;
        solid = &neigh;
        iSolid = iNeigh;
      } else {
        fluid = &neigh;
        solid = &site;
        iSolid = Neighbours::inverses[iNeigh];
      }

      LinkData& link = fluid->Links[iSolid];

      // The duct walls are always exactly half way
      link.Distance = 0.5;

      // What type of wall was crossed?
      Index type;
      int i;
      for (i = 0; i < 3; ++i) {
        // Categorise the crossing in this dimension (-1, 0,+1)
        // for (below, inside, above)
        if (solid->Position[i] < this->SquareDuct->LowerBound[i]) {
          type[i] = -1;
        } else if (solid->Position[i] > this->SquareDuct->UpperBound[i]) {
          type[i] = +1;
        } else {
          type[i] = 0;
        }
      }
      // Figure out what this means
      if (type[this->SquareDuct->OpenAxis] == 0) {
        // Def not an iolet
        link.Type = geometry::CutType::WALL;
      } else {
        // Assume it's iolet
        if (type[this->SquareDuct->OpenAxis] < 0) {
          link.Type = geometry::CutType::INLET;
          link.IoletId = 0;
        } else {
          link.Type = geometry::CutType::OUTLET;
          link.IoletId = 0;
        }

        /*
         * The wall normal at a link intersection will be normal to one
         * of the duct sides most of the times. For link intersections
         * happening exactly at the duct edge, we take the wall normal
         * to be the same as the direction of the intersecting link.
         */
        link.WallNormalAtWallCut = {0.0, 0.0, 0.0};
        for (i = 0; i < 3; ++i) {
          if (i == this->SquareDuct->OpenAxis) {
            // skip
          } else {
            if (type[i] != 0) {
              // This also crosses a wall - so wall
              link.Type = geometry::CutType::WALL;
              switch (i) {
                case 0:
                  link.WallNormalAtWallCut += Vector(type[i], 0, 0);
                  break;
                case 1:
                  link.WallNormalAtWallCut += Vector(0, type[i], 0);
                  break;
                case 2:
                  link.WallNormalAtWallCut += Vector(0, 0, type[i]);
                  break;
              }
            } else {
              // Might be iolet
            }
          }
        }
        link.WallNormalAtWallCut.Normalise();
        link.DistanceInVoxels = link.Distance * Neighbours::norms[iSolid];
      }
    }
  }

  // If there's enough information available, an approximation of the wall
  // normal will be computed for this fluid site.
  this->ComputeAveragedNormal(site);
}

}  // namespace hemelb::gmytool::gmy
