
//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_UNITTESTS_REDBLOOD_CELLARMY_H
#define HEMELB_UNITTESTS_REDBLOOD_CELLARMY_H

#include <cppunit/TestFixture.h>
#include "unittests/redblood/Fixtures.h"
#include "redblood/CellArmy.h"

namespace hemelb { namespace unittests { namespace redblood {

//! Mock cell for ease of use
class FakeCell : public hemelb::redblood::CellBase {
  public:
    mutable size_t nbcalls = 0;
    using hemelb::redblood::CellBase::CellBase;
    //! Facet bending energy
    virtual PhysicalEnergy operator()() const override { return 0; }
    //! Facet bending energy
    virtual PhysicalEnergy operator()(
        std::vector<LatticeForceVector> &) const override {
      ++nbcalls;
      return 0;
    }
};


class CellArmyTests : public helpers::FourCubeBasedTestFixture {
  CPPUNIT_TEST_SUITE(CellArmyTests);
     CPPUNIT_TEST(testCell2Fluid);
     CPPUNIT_TEST(testFluid2Cell);
  CPPUNIT_TEST_SUITE_END();

  PhysicalDistance const cutoff = 5.0;
  PhysicalDistance const halo = 2.0;
  typedef lb::lattices::D3Q15 D3Q15;
  typedef lb::kernels::GuoForcingLBGK<lb::lattices::D3Q15> Kernel;

  public:
    void testCell2Fluid();
    void testFluid2Cell();
  private:
    virtual size_t CubeSize() const { return 32 + 2; }
};

void CellArmyTests::testCell2Fluid() {
  // Fixture all pairs far from one another
  auto cells = TwoPancakeSamosas<FakeCell>(cutoff);
  assert(cells.size() == 2);
  auto cast = [](decltype(cells)::const_iterator _i)
    -> std::shared_ptr<FakeCell> {
    return std::dynamic_pointer_cast<FakeCell>(*_i);
  };
  assert(cast(cells.begin())->nbcalls == 0);
  assert(cast(std::next(cells.begin()))->nbcalls == 0);

  helpers::ZeroOutFOld(latDat);
  helpers::ZeroOutForces(latDat);

  redblood::CellArmy<Kernel> army(*latDat, cells, cutoff, halo);
  army.cell2Cell.cutoff = 0.5;
  army.cell2Cell.intensity = 1.0;
  army.cell2FluidInteractions();

  CPPUNIT_ASSERT(cast(cells.begin())->nbcalls == 1);
  CPPUNIT_ASSERT(cast(std::next(cells.begin()))->nbcalls == 1);
  for(size_t i(0); i < latDat->GetLocalFluidSiteCount(); ++i)
    CPPUNIT_ASSERT(helpers::is_zero(latDat->GetSite(i).GetForce()));

  LatticePosition const n0(15 - 0.1, 15.5, 15.5);
  LatticePosition const n1(15 + 0.2, 15.5, 15.5);
  (*cells.begin())->GetVertices().front() = n0;
  (*std::next(cells.begin()))->GetVertices().front() = n1;
  army.updateDNC();
  army.cell2FluidInteractions();

  CPPUNIT_ASSERT(cast(cells.begin())->nbcalls == 2);
  CPPUNIT_ASSERT(cast(std::next(cells.begin()))->nbcalls == 2);
  CPPUNIT_ASSERT(not helpers::is_zero(latDat->GetSite(15, 15, 15).GetForce()));
}

void CellArmyTests::testFluid2Cell() {
  // Checks that positions of cells are updated. Does not check that attendant
  // DNC is.
  auto cells = TwoPancakeSamosas<FakeCell>(cutoff);
  auto const orig = TwoPancakeSamosas<FakeCell>(cutoff);
  auto const normal = Facet(
      (*cells.begin())->GetVertices(), (*cells.begin())->GetFacets()[0]
  ).normal();

  LatticePosition gradient;
  Dimensionless non_neg_pop;
  std::function<Dimensionless(PhysicalVelocity const &)> linear, linear_inv;
  std::tie(non_neg_pop, gradient, linear, linear_inv)
    = helpers::makeLinearProfile(CubeSize(), latDat, normal);

  redblood::CellArmy<Kernel> army(*latDat, cells, cutoff, halo);
  army.fluid2CellInteractions();

  for(size_t i(0); i < cells.size(); ++i) {
    auto const disp = (*std::next(cells.begin(), i))->GetVertices().front()
      - (*std::next(orig.begin(), i))->GetVertices().front();
    auto i_nodeA = (*std::next(cells.begin(), i))->GetVertices().begin();
    auto i_nodeB = (*std::next(orig.begin(), i))->GetVertices().begin();
    auto const i_end = (*std::next(cells.begin(), i))->GetVertices().end();
    for(; i_nodeA != i_end; ++i_nodeA, ++i_nodeB)
      CPPUNIT_ASSERT(helpers::is_zero((*i_nodeA - *i_nodeB) - disp));
  }
}


CPPUNIT_TEST_SUITE_REGISTRATION(CellArmyTests);
}}}

#endif // ONCE
