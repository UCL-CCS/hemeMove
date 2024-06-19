// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "CylinderGenerator.h"

#include <cmath>
#include <queue>

#include "InconsistentFluidnessError.h"
#include "Neighbours.h"
#include "Site.h"

#include "Debug.h"

#include "io/formats/geometry.h"

namespace hemelb::gmytool::gmy {

using namespace io::formats;

CylinderGenerator::CylinderGenerator() : GeometryGenerator() {
  this->Cylinder = new CylinderData;
}

CylinderGenerator::~CylinderGenerator() {
  delete this->Cylinder;
}
void CylinderGenerator::ComputeBounds(double bounds[6]) const {
  /*
   * Compute the approximate bounds of the cylinder
   */
  Vector& c = this->Cylinder->Centre;
  Vector& n = this->Cylinder->Axis;
  double& r = this->Cylinder->Radius;
  double& h = this->Cylinder->Length;

  for (int i = 0; i < 3; ++i) {
    double sinTheta;
    double n_i = n[i] > 0. ? n[i] : -n[i];
    double n_i2 = n_i * n_i;
    if (n_i2 > 1.) {
      sinTheta = 0.;
    } else {
      sinTheta = std::sqrt(1. - n_i2);
    }
    bounds[2 * i] = c[i] - 0.5 * h * n_i - r * sinTheta;
    bounds[2 * i + 1] = c[i] + 0.5 * h * n_i + r * sinTheta;
  }
}

bool IsInsideCylinder(CylinderData* cyl, Vector& pt) {
  Vector x = pt - cyl->Centre;
  double z = Dot(x, cyl->Axis);
  if (z > -0.5 * cyl->Length && z < 0.5 * cyl->Length) {
    double rsq = Dot(x, x) - z * z;
    if (rsq < cyl->Radius * cyl->Radius) {
      return true;
    } else {
      return false;
    }
  } else {
    return false;
  }
}

#define TOL 1e-10

struct Hit {
  double t;
  Vector pt;
  int cellId;
  bool operator<(const Hit& rhs) const {
    // Note this goes the wrong way! It's meant to.
    return (t > rhs.t);
  }
};

typedef std::priority_queue<Hit> HitList;
/**
 * Compute the intersection between the line segment and the capped cylinder.
 *
 * Site from is ALWAYS fluid
 * Site to is ALWAYS solid
 *
 * Therefore we expect exactly one intersection. Due to floating point errors,
 * we can't guarantee this. So, we search for intersections with some
 * tolerance (macro TOL) and pick the closest to fluid site (as given by the
 * parametic line coordinate t. We use a priority_queue to keep the hits in
 * order.
 */
Hit ComputeIntersection(CylinderData const* cyl,
                        std::vector<Iolet> const& iolets,
                        Site& from,
                        Site& to) {
  HitList hitList = HitList();
  /*
   * The equation of the line segment is:
   * x(t) = a + t (b - a)		t E (0, 1)
   */
  auto const& n = cyl->Axis;
  auto const& r = cyl->Radius;
  auto const& c = cyl->Centre;
  auto const& h = cyl->Length;

  {
    /*
     * The equation determining intersection with an INFINITE cylinder of
     * radius r is:
     * [(b-a)^2 - ((b-a).n)^2] t^2 + [2 a.(b - a) - 2 (a.n)((b-a).n)] t + [a^2 -
     * (a.n)^2 - r^2] = 0 So first compute the coefficients and then the
     * discriminant of the eqn
     */
    Vector a = from.Position - c;
    Vector b = to.Position - c;
    Vector b_a = b - a;

    double b_aDOTn = Dot(b_a, n);
    double aDOTn = Dot(a, n);

    double A = (b_a.GetMagnitudeSquared() - b_aDOTn * b_aDOTn);
    double B = 2. * (Dot(a, b_a) - aDOTn * b_aDOTn);
    double C = a.GetMagnitudeSquared() - aDOTn * aDOTn - r * r;

    double discriminant = B * B - 4 * A * C;

    if (discriminant < 0.) {
      // No real solutions.
    } else if (discriminant == 0) {
      // Exactly one solution, i.e. line segment just brushes the cylinder.
      // This means the line must be outside the cylinder everywhere else,
      // so we will count this as no intersection.
    } else {
      // Two real solutions. So we have two intersections between the
      // infinite line and infinite cylinder.
      // If t outside (0,1), then the intersection isn't on the line segment
      // we care about.
      std::vector<double> ts(2);
      ts[0] = (-B + std::sqrt(discriminant)) / (2 * A);
      ts[1] = (-B - std::sqrt(discriminant)) / (2 * A);
      for (std::vector<double>::iterator tIt = ts.begin(); tIt != ts.end();
           ++tIt) {
        double t = *tIt;
        if (t > 0. && t < 1. + TOL) {
          // Hit in part of line we care about.
          // Now check if it's on the finite cylinder. This requires that
          // x.n E (-h/2, h/2)
          double xDOTn = aDOTn + t * b_aDOTn;
          if (xDOTn >= -0.5 * h - TOL && xDOTn <= 0.5 * h + TOL) {
            // Real cylinder hit!
            Hit hit;
            hit.t = t;
            hit.pt = from.Position + b_a * t;
            hit.cellId = -1;
            hitList.push(hit);
          }
        }
      }
    }
  }
  /*
   * Now we want to look for intersections with the capping planes.
   */
  Vector& a = from.Position;
  Vector& b = to.Position;
  Vector b_a = b - a;
  for (auto iIt = iolets.begin(); iIt != iolets.end(); ++iIt) {
    Iolet const& iolet = *iIt;
    /*
     * Plane equation is x.p = q.p (p = plane normal, q = point on plane)
     * Line is x = a + t(b-a)
     */
    Vector const& q = iolet.Centre;
    Vector const& p = iolet.Normal;

    double t = Dot(q - a, p) / Dot(b_a, p);
    if (t > 0. && t < 1. + TOL) {
      // Intersection within the line segment. Now check within cap.
      Vector x = a + b_a * t;
      Vector x_c = x - c;
      double x_cDOTn = Dot(x_c, n);
      Vector radial = x_c - n * x_cDOTn;
      if (radial.GetMagnitudeSquared() < (r + TOL) * (r + TOL)) {
        // Within the cap
        Hit hit;
        hit.t = t;
        hit.pt = x;
        hit.cellId = iIt - iolets.begin();
        hitList.push(hit);
      }
    }
  }
  // We know there SHOULD be exactly one intersection. Take the closest.
  if (hitList.size() == 0)
    throw InconsistentFluidnessError(from, to, hitList.size());

  // Ensure that the hit lies on the link
  Hit ans = hitList.top();
  if (ans.t > 1.) {
    ans.t = 1.;
    ans.pt = b;
  }
  return ans;
}

/*
 * Given a site with known fluidness, examine the links to not-yet-visited
 * neighbouring sites. If the neighbours have unknown fluidness, set that.
 *
 * Since we wish to examine each link only once, this will set link properties
 * from neighbour => site as well as site => neighbour
 *
 */
void CylinderGenerator::ClassifySite(Site& site) {
  for (LaterNeighbourIterator neighIt = site.begin(); neighIt != site.end();
       ++neighIt) {
    Site& neigh = *neighIt;
    unsigned int iNeigh = neighIt.GetNeighbourIndex();
    int nHits;

    if (!neigh.IsFluidKnown) {
      neigh.IsFluidKnown = true;
      neigh.IsFluid = IsInsideCylinder(this->Cylinder, neigh.Position);
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

      // This uses special knowledge about the fact that we're working
      // on a cylinder. Since the first point is always fluid, and the
      // second always solid, we must get exactly one hit.
      Hit hit =
          ComputeIntersection(this->Cylinder, this->Iolets, *fluid, *solid);

      LinkData& link = fluid->Links[iSolid];

      // This is set in any solid case
      float distanceInVoxels = (hit.pt - fluid->Position).GetMagnitude();
      // The distance is in voxels but must be output as a fraction of
      // the lattice vector. Scale it.
      link.Distance = distanceInVoxels / Neighbours::norms[iSolid];

      if (hit.cellId < 0) {
        // -1 => we hit a wall
        link.Type = geometry::CutType::WALL;
      } else {
        // We hit an inlet or outlet
        Iolet const& iolet = this->Iolets[hit.cellId];
        if (iolet.IsInlet) {
          link.Type = geometry::CutType::INLET;
        } else {
          link.Type = geometry::CutType::OUTLET;
        }
        // Set the Id
        link.IoletId = iolet.Id;
      }

      // If this link intersected the wall, store the normal of the cell we hit
      // and the distance to it.
      if (link.Type == geometry::CutType::WALL) {
        this->ComputeCylinderNormalAtAPoint(link.WallNormalAtWallCut, hit.pt,
                                            this->Cylinder->Axis);
        link.DistanceInVoxels = distanceInVoxels;
      }
    }
  }

  // If there's enough information available, an approximation of the wall
  // normal will be computed for this fluid site.
  this->ComputeAveragedNormal(site);
}

void CylinderGenerator::ComputeCylinderNormalAtAPoint(
    Vector& wallNormal,
    const Vector& surfacePoint,
    const Vector& cylinderAxis) const {
  wallNormal = surfacePoint - cylinderAxis * Dot(surfacePoint, cylinderAxis);
  wallNormal.Normalise();
}

}  // namespace hemelb::gmytool::gmy
