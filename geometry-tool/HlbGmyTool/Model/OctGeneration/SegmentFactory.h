// -*- mode: c++; -*-
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#ifndef HLBGMYTOOL_OCT_SEGMENTFACTORY_H
#define HLBGMYTOOL_OCT_SEGMENTFACTORY_H
#include <utility>
#include <vector>

#include "Index.h"

namespace hemelb::gmytool::oct {

// This produces all the line segments that pass through the points of a cube
// of the supplied size (one corner at the origin, the diagonal at
// (size-1, size-1, size-1)
// User needs to add the relevant offset on
class SegmentFactory {
 public:
  // A line segment from first to second
  typedef std::pair<Index, Index> Seg;

  // size = side length of the cube
  SegmentFactory(int size);
  // Make all the line segments for a given lattice vector
  std::vector<Seg> MakeSegments(const Index& vec) const;
  // Make all the starting points for a given vector
  std::vector<Index> MakeStartPoints(const Index& vec) const;

 private:
  // A single face of the axis aligned cube
  class Face {
   public:
    Face(int size, int direction, bool positive);
    std::vector<Index> MakePoints(const Index& vec) const;

   private:
    int size;
    int direction;
    bool positive;
    Index normal;
    friend SegmentFactory;
  };

  int size;
  std::vector<Face> faces;
};
}  // namespace hemelb::gmytool::oct
#endif
