// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <catch2/catch.hpp>
#include "redblood/VertexBag.h"
#include "Traits.h"
#include "tests/helpers/ApproxVector.h"

namespace hemelb
{
  namespace tests
  {
    using namespace redblood;

    using redblood::detail::procsAffectedByPosition;
    using redblood::detail::splitVertices;

    static proc_t proc_at_pos(LatticeVector const &position) {
      return (position.x() > 0 ? 0 : 1) +
	(position.y() > 0 ? 0 : 1) * 2 +
	(position.z() > 0 ? 0 : 1) * 4;
    }

    TEST_CASE("VertexBagTests", "[redblood]") {
      using Stencil = Traits<>::Stencil;

      SECTION("testConstruction") {
	auto const cell = std::make_shared<Cell>(icoSphere());
	VertexBag bag(cell);
	REQUIRE(bag.GetTag() == cell->GetTag());
	REQUIRE(size_t(0) == bag.GetVertices().size());
      }

      SECTION("testProcsAffectedBySingleProc") {
	LatticePosition const X(10, 10, 10);
	auto const high_quadrant = procsAffectedByPosition < Stencil
							     > (proc_at_pos, X);
	REQUIRE(size_t(1) == high_quadrant.size());
	REQUIRE(size_t(1) == high_quadrant.count(0));
	auto const low_quadrant = procsAffectedByPosition < Stencil
							    > (proc_at_pos, -X);
	REQUIRE(size_t(1) == low_quadrant.size());
	REQUIRE(size_t(1) == low_quadrant.count(7));

	auto const avoid = procsAffectedByPosition < Stencil
						     > (proc_at_pos, -X, 7);
	REQUIRE(size_t(0) == avoid.size());
      }

      SECTION("testProcsAffectedByMultiProcs") {
	LatticePosition const X(10, 0, 10);
	auto const two_quadrants = procsAffectedByPosition < Stencil
							     > (proc_at_pos, X);
	REQUIRE(size_t(2) == two_quadrants.size());
	REQUIRE(size_t(1) == two_quadrants.count(0));
	REQUIRE(size_t(1) == two_quadrants.count(2));

	LatticePosition const Y(0, 0, 10);
	auto const four_quadrants = procsAffectedByPosition < Stencil
							      > (proc_at_pos, Y);
	REQUIRE(size_t(4) == four_quadrants.size());
	for (proc_t i(0); i < proc_t(4); ++i) {
	  REQUIRE(size_t(1) == four_quadrants.count(i));
	}

	LatticePosition const Z(0, 0, 0);
	auto const all_quadrants = procsAffectedByPosition < Stencil
							     > (proc_at_pos, Z);
	REQUIRE(size_t(8) == all_quadrants.size());
	for (proc_t i(0); i < proc_t(8); ++i) {
	  REQUIRE(size_t(1) == all_quadrants.count(i));
	}
      }

      SECTION("testSplittingCellSingleProc") {
	auto const cell = std::make_shared<Cell>(icoSphere());
	*cell *= 5e0;
	*cell += LatticePosition(50, 50, 50) - cell->GetBarycentre();
	auto splits = splitVertices < Stencil > (proc_at_pos, cell);
	REQUIRE(size_t(1) == splits.size());
	REQUIRE(size_t(1) == splits.count(0));
	auto const &expected = cell->GetVertices();
	auto const &actual = splits[0]->GetVertices();
	REQUIRE(expected.size() == actual.size());
	for (size_t i(0); i < expected.size(); ++i) {
	  REQUIRE(ApproxV(expected[i]).Margin(1e-8) == actual[i]);
	}
      }

      SECTION("testSplittingCellMultiProcs") {
	auto const cell = std::make_shared<Cell>(icoSphere());
	*cell *= 5e0;
	*cell += LatticePosition(1, 0, 1) * 50e0 - cell->GetBarycentre();
	auto splits = splitVertices < Stencil > (proc_at_pos, cell);
	REQUIRE(size_t(2) == splits.size());
	REQUIRE(size_t(1) == splits.count(0));
	REQUIRE(size_t(1) == splits.count(2));
	// Four points are in common since their y coordinates is zero
	REQUIRE(cell->GetVertices().size() + 4 ==
		splits[0]->GetVertices().size() + splits[2]->GetVertices().size());
      }


    }

  }
}
