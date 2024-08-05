// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <catch2/catch.hpp>

#include "redblood/buffer/Columns.h"
#include "tests/redblood/Fixtures.h"
#include "tests/helpers/ApproxVector.h"

namespace hemelb::tests
{
    using namespace redblood;

    TEST_CASE("ColumnsTests", "[redblood]") {

        auto cylinder = std::make_shared<Cylinder>();
        cylinder->normal = LatticePosition(0, 0, 1);
        cylinder->origin = LatticePosition(5, 5, 0);
        cylinder->radius = 10;
        LatticePosition const colAxis{1, 0, 0};
        LatticePosition const cellAxis{1, 1, 1};
        LatticeDistance const sep = 1;

        auto is_in_cylinder = [&cylinder](LatticePosition const &a, MeshData::Vertices const &verts) -> bool {
            LatticePosition const barycentre = redblood::barycentre(verts);
            LatticePosition const n0 = cylinder->normal.GetNormalised();
            if (Dot(a, cylinder->normal) < -1e-8) {
                return false;
            }
            for (auto const &v : verts) {
                LatticePosition const x = a + v - barycentre - cylinder->origin;
                if (Cross(x, n0).GetMagnitude() >= cylinder->radius) {
                    return false;
                }
            }
            return true;
        };

        auto approx = Approx(0.0).margin(1e-8);
        auto zero = approx(0.0);
        auto one = approx(1.0);

        SECTION("testIdentityRotationMatrix") {
            auto const r = rotationMatrix(LatticePosition(1, 0, 0), LatticePosition(1, 0, 0));
            for (size_t i(0); i < 3; ++i)
                for (size_t j(0); j < 3; ++j)
                    REQUIRE(r[i][j] == (i == j ? one : zero));
        }

        SECTION("testRotationMatrix90Degrees") {
            auto const r = rotationMatrix(LatticePosition(1, 0, 0), LatticePosition(0, 1, 0));
            // a0 maps to b0
            REQUIRE(zero == (r * LatticePosition(1, 0, 0)).x());
            REQUIRE(one == (r * LatticePosition(1, 0, 0)).y());
            REQUIRE(zero == (r * LatticePosition(1, 0, 0)).z());
            // a0.Cross(b0) maps to itsel   f
            REQUIRE(zero == (r * LatticePosition(0, 0, 1)).x());
            REQUIRE(zero == (r * LatticePosition(0, 0, 1)).y());
            REQUIRE(one == (r * LatticePosition(0, 0, 1)).z());
            // a0.Cross(a0.Cross(b0)) maps to b0.Cross(a0.Cross(b0))
            REQUIRE(one == (r * LatticePosition(0, -1, 0)).x());
            REQUIRE(zero == (r * LatticePosition(0, -1, 0)).y());
            REQUIRE(zero == (r * LatticePosition(0, -1, 0)).z());
        }

        SECTION("testRotationMatrix45Degrees") {
            auto const r = rotationMatrix(LatticePosition(1, 0, 0), LatticePosition(1, 1, 0));
            LatticeDistance const sqrt2 = std::sqrt(2e0);
            // a0 maps to b0
            REQUIRE(one == (r * LatticePosition(sqrt2, 0, 0)).x());
            REQUIRE(one == (r * LatticePosition(sqrt2, 0, 0)).y());
            REQUIRE(zero == (r * LatticePosition(sqrt2, 0, 0)).z());
            // a0.Cross(b0) maps to itself
            REQUIRE(zero == (r * LatticePosition(0, 0, 1)).x());
            REQUIRE(zero == (r * LatticePosition(0, 0, 1)).y());
            REQUIRE(one == (r * LatticePosition(0, 0, 1)).z());
            // a0.Cross(a0.Cross(b0)) maps to b0.Cross(a0.Cross(b0))
            REQUIRE(zero == (r * LatticePosition(1, 1, 0)).x());
            REQUIRE(approx(sqrt2) == (r * LatticePosition(1, 1, 0)).y());
            REQUIRE(zero == (r * LatticePosition(1, 1, 0)).z());
        }

        SECTION("testMaxExtensions") {
            using buffer::detail::maxExtension;
            auto const verts = tetrahedron().GetVertices();
            LatticeDistance const s2 = std::sqrt(2e0);
            REQUIRE(approx(s2) == maxExtension(verts, LatticePosition(1, 1, 0)));
            REQUIRE(approx(s2) == maxExtension(verts, LatticePosition(-1, -1, 0)));
            REQUIRE(one == maxExtension(verts, LatticePosition(1, 0, 0)));
            REQUIRE(one == maxExtension(verts, LatticePosition(0, 0, 1)));
            REQUIRE(one == maxExtension(verts, LatticePosition(0, 0, 0.5)));
        }

        SECTION("testIterator") {
            // create deformed mesh
            MeshData::Vertices verts = tetrahedron().GetVertices();
            verts.front() += LatticePosition(0.25, 0.25, 0.25);
            buffer::ColumnPositionIterator iterator(cylinder, verts, cellAxis, colAxis, sep);

            // Find first vector: should be collinear with column
            std::vector<LatticePosition> positions;
            positions.push_back(*iterator);
            REQUIRE(is_in_cylinder(positions[0], verts));
            ++iterator;
            positions.push_back(*iterator);
            REQUIRE(is_in_cylinder(*iterator, verts));
            const auto a0 = positions[1] - positions[0];
            REQUIRE(zero == Cross(a0, colAxis).GetMagnitude());
            // going two back should be outside cylinder
            // Goin one back might not since checking size of mesh is only approximate
            REQUIRE(not is_in_cylinder(positions[0] - a0 * 2.0, verts));

            // Now go forward until we change column
            do {
                ++iterator;
                positions.push_back(*iterator);
                REQUIRE(is_in_cylinder(positions.back(), verts));
            } while ( Cross(positions.back() - positions[0], a0).GetMagnitude() < 1e-8);

            // a1 goes from last position in one column to first position in other
            auto const a1 = positions.back() - positions[positions.size() - 2];
            REQUIRE(Cross(a1, a0).GetMagnitude() > 1e-8);
            REQUIRE(Cross(Cross(a1, a0), cylinder->normal).GetMagnitude() < 1e-8);

            // Now look for first item outside current plane
            do {
                ++iterator;
                positions.push_back(*iterator);
                REQUIRE(is_in_cylinder(positions.back(), verts));
            } while ( Dot(positions.back() - positions[0], cylinder->normal) < 1e-8);

            // Check next positions are same as previous but translated along normal
            // positions.back() == positions[i] + something * cylinder->normal.
            const auto N = positions.size();
            REQUIRE(N > 1);
            for (size_t i = 0; i < 2 * N; ++i) {
                REQUIRE(approx(positions[i].x()) == positions.back().x());
                REQUIRE(approx(positions[i].y()) == positions.back().y());
                ++iterator;
                positions.push_back(*iterator);
                REQUIRE(is_in_cylinder(positions.back(), verts));
            }
        }

      SECTION("testCellDrop") {
	//using namespace hemelb::redblood::buffer;
	auto templateCell = std::make_shared<Cell>(tetrahedron());
	templateCell->moduli.bending = 5.0;
	templateCell->moduli.volume = 3.0;
	auto const rotation = rotationMatrix(cellAxis, colAxis);

	auto check = [templateCell, &rotation, &is_in_cylinder, &approx](std::shared_ptr<CellBase const> acell) {
	  auto const cell = std::static_pointer_cast<Cell const, CellBase const>(acell);
	  auto const b0 = templateCell->GetBarycentre();
	  auto const b1 = cell->GetBarycentre();
	  auto const& vertices0 = templateCell->GetVertices();
	  auto const& vertices1 = cell->GetVertices();
	  REQUIRE(templateCell->GetNumberOfNodes() == cell->GetNumberOfNodes());
	  for(site_t i(0); i < cell->GetNumberOfNodes(); ++i) {
	    // Cells are rotated with respect to original input by a given value
	    // They are translated by a value that changes at each iteration
	    // The translation is not checked here.
	    auto const expected = rotation * (vertices0[i] - b0);
	    auto const actual = vertices1[i] - b1;
	    REQUIRE(ApproxV(expected) == actual);
	    REQUIRE(is_in_cylinder(b1, cell->GetVertices()));
	  }
	  REQUIRE(templateCell->GetTemplateMesh().isSameData(cell->GetTemplateMesh()));
	  REQUIRE(approx(templateCell->moduli.bending) == cell->moduli.bending);
	  REQUIRE(approx(templateCell->moduli.volume) == cell->moduli.volume);
	  REQUIRE(approx(templateCell->moduli.surface) == cell->moduli.surface);
	};

	buffer::ColumnCellDrop dropCell(cylinder, templateCell, cellAxis, colAxis, sep);
	for (size_t i(0); i < 50; ++i) {
	  check(dropCell());
	}
      }
    }

}
