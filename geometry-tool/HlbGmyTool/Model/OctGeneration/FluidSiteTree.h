// -*- mode: c++; -*-
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#ifndef HLBGMYTOOL_OCT_FLUIDSITETREE_H
#define HLBGMYTOOL_OCT_FLUIDSITETREE_H

#include <array>
#include <memory>

#include "Oct.h"
#include "util/Vector3D.h"

namespace hemelb::gmytool::oct {

enum class Intersection { None = 0, Wall = 1, Inlet = 2, Outlet = 3 };

// This holds the data for a single link from a fluid site
struct Link {
  inline Link()
      : type(Intersection::None),
        dist(std::numeric_limits<float>::infinity()),
        id(-1) {}
  Intersection type;
  float dist;
  int id;
};

inline std::ostream& operator<<(std::ostream& os, const Link& lnk) {
  os << '(' << static_cast<int>(lnk.type) << ',' << lnk.dist << ',' << lnk.id
     << ')';
  return os;
}

template <size_t N>
std::ostream& operator<<(std::ostream& os, const std::array<Link, N>& lnk_ar) {
  auto iter = lnk_ar.begin();
  os << *iter;

  for (; iter != lnk_ar.end(); ++iter)
    os << "," << *iter;
  return os;
}

// Single precision vector
using SVector = hemelb::util::Vector3D<float>;

// This will hold the data for a single fluid site ready to be written
struct FluidSite {
  inline FluidSite() : links(), normal(), has_normal(false) {}

  std::array<Link, 26> links;
  SVector normal;
  bool has_normal;
};

typedef std::shared_ptr<FluidSite> FluidSitePtr;
struct FluidData {
  inline FluidData() : count(0), leaf() {}

  unsigned count;
  FluidSitePtr leaf;
};
typedef Octree<FluidData> FluidTree;

}  // namespace hemelb::gmytool::oct
#endif  // HLBGMYTOOL_OCT_FLUIDSITETREE_H
