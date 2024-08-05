// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <catch2/catch.hpp>

#include "redblood/GridAndCell.h"
#include "redblood/Mesh.h"

#include "tests/helpers/ApproxVector.h"
#include "tests/helpers/LatticeDataAccess.h"
#include "tests/redblood/Fixtures.h"

namespace hemelb::tests
{
    using namespace redblood;

    namespace
    {
      // Creates a mesh that is a single planar triangle
      // It still has two facets so that we can run fake forces on it.
      Mesh triangleMesh()
      {
	// Rotate samosa so it is in xy plane
	Mesh result = pancakeSamosa(0);
	result.GetData()->vertices[1] = LatticePosition(0, 1, 0);
	result.GetData()->vertices[2] = LatticePosition(std::sqrt(3.0) / 2.0, 0.5, 0);
	return result;
      }
    }

    // Because Catch doesn't support template test cases + fixture,
    // only test cases over a templated fixture.
    template <typename STENCIL>
    class CellForceSpreadTestsFixture : public SquareDuctTetrahedronFixture {
    public:
      std::vector<LatticeForceVector> forces;
      LatticePosition center;

      CellForceSpreadTestsFixture() : SquareDuctTetrahedronFixture{triangleMesh(), 3}, center{double(cubeSizeWithHalo / 2)}
      {
      }

      LatticeForceVector force_at_center(LatticePosition const &position)
      {
        mesh += position - mesh.GetBarycentre();
        helpers::ZeroOutForces(*latDat);
        details::spreadForce2Grid<details::SpreadForces, STENCIL>(std::shared_ptr<CellBase>(&mesh,
                                                                                            [](CellBase*)
                                                                                            {}),
                                                                  details::SpreadForces(forces,
                                                                                        *latDat));
	// TODO #759: is this truncation OK?
        return latDat->GetSite(LatticeVector{center}).GetForce();
      }
    };

    TEMPLATE_TEST_CASE_METHOD(CellForceSpreadTestsFixture,
				   "CellForceSpreadTests",
				   "[redblood]",
				   stencil::FourPoint, stencil::CosineApprox, stencil::ThreePoint, stencil::TwoPoint) {
      using STENCIL = TestType;

      LatticePosition direction;
      auto setupForces = [this, &direction]() {
	direction = LatticePosition(1, 2, 3);
	LatticePosition intensity = LatticePosition(3, 2, 1).Normalise();
	this->forces.resize(this->mesh.GetNumberOfNodes());
	auto i_vertex = this->mesh.GetVertices().begin();
	auto const i_end = this->mesh.GetVertices().end();
	auto i_force = this->forces.begin();
      
	for (; i_vertex != i_end; ++i_vertex, ++i_force) {
	  *i_force = direction * Dot(*i_vertex, intensity);
	}
      };
      setupForces();

      SECTION("IsZeroFarFromMembrane") {
        // Very far away from pancake samosa
        auto const border = Dimensionless(STENCIL::GetRange()) / 2e0;
	
        LatticeForceVector const faraway = this->force_at_center(this->center
								 + LatticePosition(0, 0, border + 1e0));

      	auto zero = ApproxVector<LatticeForceVector>{0}.Margin(1e-8);
        REQUIRE(faraway == zero);

        auto const justTooFar = this->force_at_center(this->center + LatticePosition(0, 0, border));
        REQUIRE(justTooFar == zero);

        auto const justInside = this->force_at_center(this->center + LatticePosition(0, 0, border - 1e-4));
        REQUIRE(justInside != zero);
        REQUIRE(justInside[0] > 1e-8);
      	REQUIRE(justInside[0] < 1e-4);
      }

      SECTION("testIsSymmetric") {
        size_t const N = 10;

        for (size_t i = 0; i < N; ++i)
        {
          auto const range = Dimensionless(STENCIL::GetRange() * 0.5);
          auto const disp = range * (1e0 - Dimensionless(i) / Dimensionless(N));
          LatticePosition const displacement(0, 0, disp);
          LatticeForceVector const left = this->force_at_center(this->center + displacement);
          LatticeForceVector const right = this->force_at_center(this->center - displacement);
          REQUIRE(left == ApproxV(right));
        }
      }

      SECTION("testIsIncreasing") {
        size_t const N(10);
        LatticeForceVector last(0, 0, 0);

        for (size_t i(0); i < N; ++i) {
	  auto const range = Dimensionless(STENCIL::GetRange() * 0.5);
	  auto const d = range * (1e0 - Dimensionless(i + 1) / Dimensionless(N));
	  LatticePosition const displacement(0, 0, d);
	  LatticeForceVector const current = this->force_at_center(this->center + displacement);
	  REQUIRE(current[0] > last[0]);
	  REQUIRE(current[1] > last[1]);
	  REQUIRE(current[2] > last[2]);
	  last = current;
	}
      }

      SECTION("testIsLinear") {
        size_t const N(5);
        this->mesh = Cell(refine(MeshData { this->mesh.GetVertices(), this->mesh.GetFacets() }, 4));
        setupForces();
        // x0, x1 should be further than 2 from the edges
        // Only linear if samosa appears as infinite plane
        // with sufficiently dense vertices
	auto& center = this->center;
        LatticePosition const x0(center[0], center[1] - 0.5, center[2] - 0.1);
        LatticePosition const x1(center[0], center[1] + 0.5, center[2] - 0.1);
        LatticeForceVector const v0(this->force_at_center(x0)),
	  v1(this->force_at_center(x1));

        LatticeForceVector const a( (v1 - v0) / (Dot(direction, x1) - Dot(direction, x0)));

        auto const reltol = std::vector<double> { 1e-3, 1e-4, 1e-5 }.at(STENCIL::GetRange() - 2);
	auto approx = Approx(0.0).epsilon(reltol).margin(1e-8);
        for (size_t i(0); i < N; ++i)
        {
          LatticePosition const x = (x1 - x0) * (Dimensionless(i + 1) / Dimensionless(N + 2)) + x0;
          LatticeForceVector const expected(a * (Dot(direction, x) - Dot(direction, x0)) + v0);
          auto const actual = this->force_at_center(x);
          REQUIRE(approx(expected.x()) == actual.x());
	  REQUIRE(approx(expected.y()) == actual.y());
	  REQUIRE(approx(expected.z()) == actual.z());
        }
      }


    }
}
