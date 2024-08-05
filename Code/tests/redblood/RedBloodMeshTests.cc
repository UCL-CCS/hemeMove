// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <sstream>
#include <catch2/catch.hpp>
#include "resources/Resource.h"
#include "redblood/Mesh.h"
#include "redblood/MeshIO.h"

#include "tests/redblood/Fixtures.h"
#include "tests/helpers/ApproxVector.h"

namespace hemelb
{
  namespace tests
  {
    using namespace hemelb::redblood;
    // Tests geometry properties
    TEST_CASE_METHOD(BasisFixture, "RedBloodMeshTests") {
      redblood::KruegerMeshIO io = {};
      mesh.vertices.at(0) = LatticePosition(0.1, 0.1, 0.1);

      auto checkMeshOrientation = [](MeshData const &mesh) {
	MeshData copy(mesh);
	orientFacets(copy);
	for (size_t i(0); i < mesh.facets.size(); ++i) {
	  REQUIRE(mesh.facets[i][0] == copy.facets[i][0]);
	  REQUIRE(mesh.facets[i][1] == copy.facets[i][1]);
	  REQUIRE(mesh.facets[i][2] == copy.facets[i][2]);
	}
      };

      SECTION("testOrientation") {
	checkMeshOrientation(mesh);

	// change facet order of copy for next test
	MeshData copy(mesh);
	for (auto &facet : copy.facets) {
	  std::swap(facet[0], facet[2]);
	}
	// re-orient
	orientFacets(copy);
	checkMeshOrientation(copy);
      }

      SECTION("testStandardTestMeshOrientation") {
	//const auto io = redblood::KruegerMeshIO{};
	checkMeshOrientation(*io.readFile(resources::Resource("red_blood_cell.txt").Path(), true));
	checkMeshOrientation(*tetrahedron(2).GetData());
	checkMeshOrientation(*icoSphere(5).GetData());
      }

      SECTION("testVolume") {
	LatticePosition a(0.1, 0.1, 0.1), b(1, 0, 0), c(0, 1, 0), d(0, 0, 1);
	double const expected = std::abs( Dot(Cross(b - a, c - a), d - a)) / 6.0;

	for (size_t facet(0); facet != 4; ++facet) {
	  REQUIRE(Approx(expected).margin(1e-8) == volume(mesh));
	  // swap vertices so orientation stays the same, check volume still
	  // correct. This makes sure that volume does not depend on the order of
	  // the vertices. It does depend on facet orientation. However, so do
	  // the forces on the particle.
	  size_t const v0(mesh.facets[facet][0]), v1(mesh.facets[facet][1]),
	    v2(mesh.facets[facet][2]);
	  mesh.facets[facet][0] = v1;
	  mesh.facets[facet][1] = v2;
	  mesh.facets[facet][2] = v0;
	  REQUIRE(Approx(expected).margin(1e-8) == volume(mesh));

	  mesh.facets[facet][0] = v2;
	  mesh.facets[facet][1] = v0;
	  mesh.facets[facet][2] = v1;
	  REQUIRE(Approx(expected).margin(1e-8) == volume(mesh));

	  // put back to first order
	  mesh.facets[facet][0] = v0;
	  mesh.facets[facet][1] = v1;
	  mesh.facets[facet][2] = v2;
	}
      }

      SECTION("testBarycentre") {
	LatticePosition a(0.1, 0.1, 0.1), b(1, 0, 0), c(0, 1, 0), d(0, 0, 1);
	LatticePosition const expected = (a + b + c + d) * 0.25;
	REQUIRE(std::abs(barycentre(mesh)[0] - expected[0]) < 1e-8);
	REQUIRE(std::abs(barycentre(mesh)[1] - expected[1]) < 1e-8);
	REQUIRE(std::abs(barycentre(mesh)[2] - expected[2]) < 1e-8);
      }

      SECTION("testScaling") {
	Dimensionless const scale = 2.5;
	Mesh original(mesh);
	Mesh scaled(mesh);
	scaled *= scale;

	REQUIRE(ApproxV(original.GetBarycentre()) == scaled.GetBarycentre());
	LatticePosition const first = (*original.GetVertices().begin()
				       - original.GetBarycentre()) * scale + original.GetBarycentre();
	LatticePosition const second = (* (++original.GetVertices().begin())
					- original.GetBarycentre()) * scale + original.GetBarycentre();
	REQUIRE(ApproxV(first) == *scaled.GetVertices().begin());
	REQUIRE(ApproxV(second) ==  *(++scaled.GetVertices().begin()));
      }

      SECTION("testTranslation") {
	LatticePosition const offset(1, 2, 3);
	Mesh original(mesh);
	Mesh trans(mesh);
	trans += offset;

	REQUIRE(ApproxV(original.GetBarycentre() + offset) ==
		trans.GetBarycentre());
	LatticePosition const first = *original.GetVertices().begin() + offset;
	LatticePosition const second = * (++original.GetVertices().begin()) + offset;
	REQUIRE(ApproxV(first) == *trans.GetVertices().begin());
	REQUIRE(ApproxV(second) == *(++trans.GetVertices().begin()));
      }

      SECTION("testResourceVolume") {
	auto const path = resources::Resource("red_blood_cell.txt").Path();
	REQUIRE(volume(*io.readFile(path, true)) > 0e0);
      }
    }
  }
}

