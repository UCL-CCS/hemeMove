//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_UNITTESTS_REDBLOOD_PARTICLEIMLP_H
#define HEMELB_UNITTESTS_REDBLOOD_PARTICLEIMLP_H

#include <cppunit/TestFixture.h>
#include "resources/Resource.h"
#include "redblood/ParticleImpl.cc"
#include "fixtures.h"

namespace hemelb { namespace unittests {

// Tests functionality that is *not* part of the HemeLB API
// Checks that we know how to compute geometric properties between facets
// However, HemeLB only cares about energy and forces
class FacetTests : public TetrahedronFixture {
  CPPUNIT_TEST_SUITE(FacetTests);
  CPPUNIT_TEST(testNormal);
  CPPUNIT_TEST(testUnitNormal);
  CPPUNIT_TEST(testAngle);
  CPPUNIT_TEST(testCommonNodes);
  CPPUNIT_TEST(testOrientedAngle);
  CPPUNIT_TEST(testSingleNodes);
  CPPUNIT_TEST_SUITE_END();
public:
  void setUp() {
    TetrahedronFixture::setUp();
    main.reset(new redblood::Facet(mesh, 0));
    neighbor.reset(new redblood::Facet(mesh, 3));
  }

  void tearDown() {}

  void testNormal() {
    {
      LatticePosition const actual = normal(*main);
      LatticePosition const expected(0, 0, -1);
      CPPUNIT_ASSERT(is_zero(actual - expected));
    }
    {
      LatticePosition const actual = normal(*neighbor);
      LatticePosition const expected(1, 1, 1);
      CPPUNIT_ASSERT(is_zero(actual - expected));
    }
  }

  void testUnitNormal() {
    {
      LatticePosition const actual = unitNormal(*main);
      LatticePosition const expected(0, 0, -1);
      CPPUNIT_ASSERT(is_zero(actual - expected));
    } {
      LatticePosition const actual = unitNormal(*neighbor);
      LatticePosition const expected(1, 1, 1);
      CPPUNIT_ASSERT(is_zero(actual - expected.GetNormalised()));
    }
  }

  void testAngle() {
    Angle const actual0 = redblood::angle(*main, *neighbor);
    CPPUNIT_ASSERT(is_zero(actual0 - std::acos(-1.0 /std::sqrt(3.))));

    mesh.vertices.back()[2] = 1e0 / std::sqrt(2.0);
    Angle const actual1 = redblood::angle(*main, *neighbor);
    CPPUNIT_ASSERT(is_zero(actual1 - 3.0*PI/4.0));
    mesh.vertices.back()[1] = 1e0;
  }

  void testCommonNodes() {
    redblood::t_IndexPair const nodes \
        = redblood::commonNodes(*main, *neighbor);
    CPPUNIT_ASSERT(nodes.first == 1 and nodes.second == 2);
    LatticePosition const edge = commonEdge(*main, *neighbor);
    CPPUNIT_ASSERT(is_zero(edge - (mesh.vertices[1] - mesh.vertices[2])));
  }
  void testSingleNodes() {
    redblood::t_IndexPair const nodes 
      = redblood::singleNodes(*main, *neighbor);
    // indices into the list of vertices of each facet,
    // not into list of all vertices.
    CPPUNIT_ASSERT(nodes.first == 0 and nodes.second == 1);
  }

  void testOrientedAngle() {
    Angle const actual0 = redblood::orientedAngle(*main, *neighbor);
    CPPUNIT_ASSERT(is_zero(actual0 + std::acos(-1.0 /std::sqrt(3.))));

    // simpler angle
    mesh.vertices.back()[2] = 1e0 / std::sqrt(2.0);
    std::cout << "\n";
    Angle const actual1 = redblood::orientedAngle(*main, *neighbor);
    CPPUNIT_ASSERT(is_zero(actual1 + 3.0*PI/4.0));

    // change orientation <==> negative angle
    mesh.facets.front()[1] = 2;
    mesh.facets.front()[2] = 1;
    Angle const actual2 = redblood::orientedAngle(
        redblood::Facet(mesh, 0), redblood::Facet(mesh, 3)
    );
    CPPUNIT_ASSERT(is_zero(actual2 + PI/4.0));
    mesh.facets.front()[1] = 1;
    mesh.facets.front()[2] = 2;

    mesh.vertices.back()[1] = 1e0;
  }
protected:
  LatticePosition nodes[4];
  std::set<size_t> main_indices;
  std::set<size_t> neighbor_indices;
  boost::shared_ptr<redblood::Facet> main, neighbor;
};

class EnergyTests : public TetrahedronFixture {
    CPPUNIT_TEST_SUITE(EnergyTests);
    CPPUNIT_TEST(testBending);
    CPPUNIT_TEST(testVolume);
    CPPUNIT_TEST(testSurface);
    CPPUNIT_TEST_SUITE_END();
  public:
    void setUp() {
      TetrahedronFixture::setUp();
      original = mesh;
      forces.resize(4, LatticeForceVector(0, 0, 0));
    }

    void tearDown() { TetrahedronFixture::tearDown(); }

    void testBending() {
      // No difference between original and current mesh
      // Hence energy is zero
      PhysicalEnergy const actual0(
        redblood::facetBending(mesh, original, 0, 3, 1e0)
      );
      CPPUNIT_ASSERT(is_zero(actual0));

      // Now modify mesh and check "energy" is square of angle difference
      mesh.vertices.back()[2] = 1e0 / std::sqrt(2.0);
      PhysicalEnergy const actual1(
        redblood::facetBending(mesh, original, 0, 3, 1e0)
      );
      mesh.vertices.back()[2] = 1e0;

      PhysicalEnergy const expected(
          std::pow((PI / 4e0 - std::acos(1./std::sqrt(3.))), 2)
      );
      CPPUNIT_ASSERT(is_zero(actual1 - expected));
    }

    void testVolume() {
      // No difference between original and current mesh
      // Hence energy is zero
      PhysicalEnergy const actual0(
        redblood::volumeEnergy(mesh, original, 1e0)
      );
      CPPUNIT_ASSERT(is_zero(actual0));

      // Now modify mesh and check "energy" is square of volume diff
      mesh.vertices.back()[2] = 1e0 / std::sqrt(2.0);
      PhysicalEnergy const actual1(
        redblood::volumeEnergy(mesh, original, 2.0*volume(original))
      );

      PhysicalEnergy const deltaV(volume(mesh) - volume(original));
      CPPUNIT_ASSERT(is_zero(actual1 - deltaV * deltaV));
      mesh.vertices.back()[2] = 1e0;
    }

    void testSurface() {
      // No difference between original and current mesh
      // Hence energy is zero
      PhysicalEnergy const actual0(
        redblood::surfaceEnergy(mesh, original, 1e0)
      );
      CPPUNIT_ASSERT(is_zero(actual0));

      // Now modify mesh and check "energy" is square of volume diff
      mesh.vertices.back()[2] = 1e0 / std::sqrt(2.0);
      PhysicalEnergy const actual1(
        redblood::surfaceEnergy(mesh, original, 2.0*surface(original))
      );

      PhysicalEnergy const deltaS(surface(mesh) - surface(original));
      CPPUNIT_ASSERT(is_zero(actual1 - deltaS * deltaS));
      mesh.vertices.back()[2] = 1e0;
    }

  protected:
    redblood::MeshData original;
    std::vector<LatticeForceVector> forces;
};



struct EnergyFunctional {
  EnergyFunctional(redblood::MeshData const &_original)
    : original(_original) {};
  redblood::MeshData const &original;
};
#define HEMELB_OP_MACRO(CLASS, FUNCTION)                                     \
  struct CLASS : public EnergyFunctional {                                   \
    CLASS(redblood::MeshData const &_original)                               \
      : EnergyFunctional(_original) {}                                       \
    PhysicalEnergy operator()(redblood::MeshData const &_mesh) const {       \
      return redblood::FUNCTION(_mesh, original, 1.0);                       \
    }                                                                        \
    PhysicalEnergy operator()(redblood::MeshData const &_mesh,               \
      std::vector<LatticeForceVector> &_forces) const{                       \
      return redblood::FUNCTION(_mesh, original, 1.0, _forces);              \
    }                                                                        \
  }
HEMELB_OP_MACRO(VolumeEnergy, volumeEnergy);
HEMELB_OP_MACRO(SurfaceEnergy, surfaceEnergy);
#undef HEMELB_OP_MACRO

struct BendingEnergy : public EnergyFunctional {
  BendingEnergy(redblood::MeshData const &_original)
    : EnergyFunctional(_original) {}
  PhysicalEnergy operator()(redblood::MeshData const &_mesh) const {
    return redblood::facetBending(_mesh, original, 0, 1, 1.0)
      + redblood::facetBending(_mesh, original, 0, 2, 1.0)
      + redblood::facetBending(_mesh, original, 0, 3, 1.0)
      + redblood::facetBending(_mesh, original, 1, 2, 1.0)
      + redblood::facetBending(_mesh, original, 1, 3, 1.0)
      + redblood::facetBending(_mesh, original, 2, 3, 1.0);
  }
  PhysicalEnergy operator()(redblood::MeshData const &_mesh,
    std::vector<LatticeForceVector> &_forces) const{
    return redblood::facetBending(_mesh, original, 0, 1, 1.0, _forces)
      + redblood::facetBending(_mesh, original, 0, 2, 1.0, _forces)
      + redblood::facetBending(_mesh, original, 0, 3, 1.0, _forces)
      + redblood::facetBending(_mesh, original, 1, 2, 1.0, _forces)
      + redblood::facetBending(_mesh, original, 1, 3, 1.0, _forces)
      + redblood::facetBending(_mesh, original, 2, 3, 1.0, _forces);
  }
};


class GradientTests : public EnergyVsGradientFixture {
    CPPUNIT_TEST_SUITE(GradientTests);
    CPPUNIT_TEST(testBending);
    CPPUNIT_TEST(testSurface);
    CPPUNIT_TEST(testVolume);
    CPPUNIT_TEST_SUITE_END();
public:

    void setUp() {
      TetrahedronFixture::setUp();
      original = mesh;
      mesh.vertices[0] += LatticePosition(-0.01, 0.02342, 0.03564);
      mesh.vertices[1] += LatticePosition(0.0837, -0.012632, 0.0872935);
      mesh.vertices[2] += LatticePosition(0.02631, -0.00824223, -0.098362);
    }

    void testBending() { energyVsForces(BendingEnergy(original)); }
    void testSurface() { energyVsForces(SurfaceEnergy(original)); }
    void testVolume() { energyVsForces(VolumeEnergy(original)); }

protected:
    redblood::MeshData original;
};

// Tests certain node movement that do not result in forces/energy
class GradientKernTests : public TetrahedronFixture {
    CPPUNIT_TEST_SUITE(GradientKernTests);
    CPPUNIT_TEST(testBending);
    CPPUNIT_TEST(testVolume);
    CPPUNIT_TEST(testSurface);
    CPPUNIT_TEST_SUITE_END();

#   define HEMELB_MACRO(NAME, FUNCTION_CALL)                                  \
    void NAME(LatticePosition const &_normal, size_t _node) {                 \
      std::vector<LatticeForceVector> forces(4, LatticeForceVector(0, 0, 0)); \
      double const epsilon(1e-5);                                             \
      redblood::MeshData newmesh(mesh);                                       \
      newmesh.vertices[_node] += _normal * epsilon;                           \
      PhysicalEnergy const deltaE(redblood::FUNCTION_CALL);                   \
      CPPUNIT_ASSERT(is_zero(deltaE / epsilon));                              \
      for(size_t i(0); i < forces.size(); ++i)                                \
        CPPUNIT_ASSERT(is_zero(forces[i]));                                   \
    }

      HEMELB_MACRO(noBending, facetBending(newmesh, mesh, 0, 3, 1.0, forces));
      HEMELB_MACRO(noVolume, volumeEnergy(newmesh, mesh, 1.0, forces));
      HEMELB_MACRO(noSurface, surfaceEnergy(newmesh, mesh, 1.0, forces));
#   undef HEMELB_MACRO

    void testBending() {
      // Single nodes
      noBending(LatticeForceVector(1, 0, 0).GetNormalised(), 0);
      noBending(LatticeForceVector(0, 1, 0).GetNormalised(), 0);
      noBending(LatticeForceVector(1, 1, -2).GetNormalised(), 3);
      noBending(LatticeForceVector(1, -1, 0).GetNormalised(), 3);
      // Common nodes
      noBending(LatticeForceVector(1, -1, 0).GetNormalised(), 1);
      noBending(LatticeForceVector(1, -1, 0).GetNormalised(), 2);
    }

    void testVolume() {
      noVolume(LatticePosition(-1, 1, 0).GetNormalised(), 0);
      noVolume(LatticePosition(1, 1, -2).GetNormalised(), 0);
      noVolume(LatticePosition(0, 1, 0).GetNormalised(), 1);
      noVolume(LatticePosition(0, 0, 1).GetNormalised(), 1);
      noVolume(LatticePosition(1, 0, 0).GetNormalised(), 2);
      noVolume(LatticePosition(0, 0, 1).GetNormalised(), 2);
      noVolume(LatticePosition(1, 0, 0).GetNormalised(), 3);
      noVolume(LatticePosition(0, 1, 0).GetNormalised(), 3);
    }

    void testSurface() {
      // Get rid of some surfaces first, to simplify kernel
      redblood::MeshData const save(mesh);
      mesh.facets.resize(2);
      mesh.facets.back() = save.facets.back();

      try {
        noSurface(LatticePosition(0, 0, 1).GetNormalised(), 0);
        noSurface(LatticePosition(1, 1, 1).GetNormalised(), 3);
        noSurface(LatticePosition(-1, 1, 0).GetNormalised(), 0);
        noSurface(LatticePosition(-1, 1, 0).GetNormalised(), 3);
      } catch(...) {
        mesh = save;
        throw;
      }
      mesh = save;
    }
};

CPPUNIT_TEST_SUITE_REGISTRATION(FacetTests);
CPPUNIT_TEST_SUITE_REGISTRATION(EnergyTests);
CPPUNIT_TEST_SUITE_REGISTRATION(GradientTests);
CPPUNIT_TEST_SUITE_REGISTRATION(GradientKernTests);
}}

#endif // ONCE
