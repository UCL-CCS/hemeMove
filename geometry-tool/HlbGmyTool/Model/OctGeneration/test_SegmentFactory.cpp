// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <catch2/catch.hpp>
#include "Neighbours.h"
#include "SegmentFactory.h"

namespace hemelb::gmytool::oct {

class SegmentFactoryTests {
  using Seg = SegmentFactory::Seg;

 public:
  void SinglePoints() {
    // Check that for a single point we create everything correctly
    SegmentFactory sf(1);
    auto&& dirs = Neighbours::Displacements;

    for (auto dir : dirs) {
      auto points = sf.MakeStartPoints(dir);
      // Only one
      REQUIRE(points.size() == 1);
      // Should start in the opposite direction
      REQUIRE(points[0] == -dir);
    }
  }
  void SingleSegs() {
    // Check that for a single point we create everything correctly
    SegmentFactory sf(1);
    auto&& dirs = Neighbours::Displacements;

    for (auto dir : dirs) {
      auto segs = sf.MakeSegments(dir);
      // Only one
      REQUIRE(segs.size() == 1);
      // Should start in the opposite direction
      REQUIRE(segs[0].first == -dir);
      REQUIRE(segs[0].second == dir);
    }
  }

  void Four() {
    SegmentFactory sf(4);
    auto&& dirs = Neighbours::Displacements;
    Index box_min(0);
    Index box_max(3);

    for (auto dir : dirs) {
      auto segs = sf.MakeSegments(dir);

      // can be 1,2,3 based on direction type
      int mag = dir.GetMagnitudeSquared();
      // Use the type to figure out how many there should be
      switch (mag) {
        case 1:
          //<100> type
          REQUIRE(segs.size() == 16);
          break;
        case 2:
          // <110> type
          REQUIRE(segs.size() == (16 + 12));
          break;
        case 3:
          // <111> type
          REQUIRE(segs.size() == (16 + 12 + 9));
          break;
        default:
          FAIL("Direction of unknown type!");
          break;
      }

      // Check that all points between first and second are in the box
      for (auto seg : segs) {
        auto start = seg.first;
        auto end = seg.second;
        for (auto pt = start + dir; !(pt == end); pt += dir) {
          REQUIRE(pt.IsInRange(box_min, box_max));
        }
      }
    }
  }
};

METHOD_AS_TEST_CASE(SegmentFactoryTests::SinglePoints,
                    "SinglePoints",
                    "[segment]");
METHOD_AS_TEST_CASE(SegmentFactoryTests::SingleSegs, "SingleSegs", "[segment]");
METHOD_AS_TEST_CASE(SegmentFactoryTests::Four, "Four", "[segment]");
}  // namespace hemelb::gmytool::oct
