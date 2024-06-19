// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#include "SegmentFactory.h"
#include <cassert>

namespace hemelb::gmytool::oct {

SegmentFactory::Face::Face(int s, int d, bool p)
    : size(s), direction(d), positive(p), normal() {
  assert(direction >= 0 && direction <= 2);
  if (positive)
    normal[direction] = 1;
  else
    normal[direction] = -1;
}

std::vector<Index> SegmentFactory::Face::MakePoints(const Index& vec) const {
  assert(Index::Dot(vec, normal) > 0);
  Index lo = Index::Zero() - vec;
  Index hi = lo + Index{size};

  lo[direction] = positive ? -1 : size;
  hi[direction] = lo[direction] + 1;

  std::vector<Index> pset;
  pset.reserve(size * size);
  for (int i = lo.x; i < hi.x; ++i)
    for (int j = lo.y; j < hi.y; ++j)
      for (int k = lo.z; k < hi.z; ++k) {
        pset.emplace_back(i, j, k);
      }
  return pset;
}

SegmentFactory::SegmentFactory(int s) : size(s) {
  faces = {Face(size, 0, true),  Face(size, 0, false), Face(size, 1, true),
           Face(size, 1, false), Face(size, 2, true),  Face(size, 2, false)};
}

std::vector<Index> SegmentFactory::MakeStartPoints(const Index& vec) const {
  std::vector<Index> start_points;

  for (auto f : faces) {
    if (Index::Dot(f.normal, vec) > 0) {
      auto face_pts = f.MakePoints(vec);
      // Pre-reserve space so we don't invalidate the iterators
      start_points.reserve(start_points.size() + face_pts.size());
      auto begin = start_points.begin();
      auto end = start_points.end();

      for (auto pt : face_pts) {
        if (std::find(begin, end, pt) == end)
          start_points.push_back(pt);
      }
    }
  }
  return start_points;
}

auto SegmentFactory::MakeSegments(const Index& vec) const -> std::vector<Seg> {
  std::vector<Seg> segments;

  auto in_box = [&](const Index& pt) {
    return (pt.x >= 0 && pt.y >= 0 && pt.z >= 0) &&
           (pt.x < size && pt.y < size && pt.z < size);
  };

  auto start_points = MakeStartPoints(vec);
  segments.reserve(start_points.size());
  for (auto start : start_points) {
    // Going to step along the direction
    // Always going to need at least 2
    auto end = start + vec * 2;
    while (in_box(end))
      end += vec;

    segments.emplace_back(start, end);
  }
  return segments;
}

}  // namespace hemelb::gmytool::oct
