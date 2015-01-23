
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

namespace hemelb
{
  namespace unittests
  {
    namespace redblood
    {

      //! Mock cell for ease of use
      class FakeCell : public hemelb::redblood::CellBase
      {
        public:
          mutable size_t nbcalls = 0;
          using hemelb::redblood::CellBase::CellBase;
          //! Facet bending energy
          virtual PhysicalEnergy operator()() const override
          {
            return 0;
          }
          //! Facet bending energy
          virtual PhysicalEnergy operator()(
            std::vector<LatticeForceVector> &) const override
          {
            ++nbcalls;
            return 0;
          }
      };


      class CellArmyTests : public helpers::FourCubeBasedTestFixture
      {
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
          virtual size_t CubeSize() const
          {
            return 32 + 2;
          }
      };

      void CellArmyTests::testCell2Fluid()
      {
        // Fixture all pairs far from one another
        auto cells = TwoPancakeSamosas<FakeCell>(cutoff);
        assert(cells.size() == 2);
        assert(std::dynamic_pointer_cast<FakeCell>(cells.front())->nbcalls == 0);
        assert(std::dynamic_pointer_cast<FakeCell>(cells.back())->nbcalls == 0);

        helpers::ZeroOutFOld(latDat);
        helpers::ZeroOutForces(latDat);

        redblood::CellArmy<Kernel> army(*latDat, cells, cutoff, halo);
        army.cell2Cell.cutoff = 0.5;
        army.cell2Cell.intensity = 1.0;
        army.cell2FluidInteractions();

        CPPUNIT_ASSERT(
          std::dynamic_pointer_cast<FakeCell>(cells.front())->nbcalls == 1);
        CPPUNIT_ASSERT(
          std::dynamic_pointer_cast<FakeCell>(cells.back())->nbcalls == 1);

        for(size_t i(0); i < latDat->GetLocalFluidSiteCount(); ++i)
        {
          CPPUNIT_ASSERT(helpers::is_zero(latDat->GetSite(i).GetForce()));
        }

        LatticePosition const n0(15 - 0.1, 15.5, 15.5);
        LatticePosition const n1(15 + 0.2, 15.5, 15.5);
        cells.front()->GetVertices().front() = n0;
        cells.back()->GetVertices().front() = n1;
        army.updateDNC();
        army.cell2FluidInteractions();

        CPPUNIT_ASSERT(
          std::dynamic_pointer_cast<FakeCell>(cells.front())->nbcalls == 2);
        CPPUNIT_ASSERT(
          std::dynamic_pointer_cast<FakeCell>(cells.back())->nbcalls == 2);
        CPPUNIT_ASSERT(not helpers::is_zero(latDat->GetSite(15, 15, 15).GetForce()));
      }

      void CellArmyTests::testFluid2Cell()
      {
        // Checks that positions of cells are updated. Does not check that attendant
        // DNC is.
        auto cells = TwoPancakeSamosas<FakeCell>(cutoff);
        auto const orig = TwoPancakeSamosas<FakeCell>(cutoff);
        auto const normal = Facet(
                              cells[0]->GetVertices(), cells[0]->GetFacets()[0]
                            ).normal();

        LatticePosition gradient;
        Dimensionless non_neg_pop;
        std::function<Dimensionless(PhysicalVelocity const &)> linear, linear_inv;
        std::tie(non_neg_pop, gradient, linear, linear_inv)
          = helpers::makeLinearProfile(CubeSize(), latDat, normal);

        redblood::CellArmy<Kernel> army(*latDat, cells, cutoff, halo);
        army.fluid2CellInteractions();

        for(size_t i(0); i < cells.size(); ++i)
        {
          auto const disp = cells[i]->GetVertices().front()
                            - orig[i]->GetVertices().front();
          auto i_nodeA = cells[i]->GetVertices().begin();
          auto i_nodeB = orig[i]->GetVertices().begin();
          auto const i_end = cells[i]->GetVertices().end();

          for(; i_nodeA != i_end; ++i_nodeA, ++i_nodeB)
          {
            CPPUNIT_ASSERT(helpers::is_zero((*i_nodeA - *i_nodeB) - disp));
          }
        }
      }


      CPPUNIT_TEST_SUITE_REGISTRATION(CellArmyTests);
    }
  }
}

#endif // ONCE
