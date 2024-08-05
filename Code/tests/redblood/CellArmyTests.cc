// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <catch2/catch.hpp>

#include "redblood/CellArmy.h"
#include "redblood/Facet.h"
#include "Traits.h"

#include "tests/redblood/Fixtures.h"
#include "tests/helpers/ApproxVector.h"
#include "tests/helpers/LatticeDataAccess.h"

namespace hemelb::tests
{
    //! Mock cell for ease of use
    class FakeCell : public redblood::Cell
    {
    public:
        mutable size_t nbcalls = 0;
        using redblood::Cell::Cell;
        //! Facet bending energy
        LatticeEnergy operator()() const override
        {
            return 0;
        }
        //! Facet bending energy
        LatticeEnergy operator()(std::vector<LatticeForceVector> &) const override
        {
            ++nbcalls;
            return 0;
        }
    };

    TEST_CASE_METHOD(helpers::FourCubeBasedTestFixture<32>, "CellArmyTests", "[redblood]") {
      using namespace redblood;

      LatticeDistance const cutoff = 5.0;
        using Traits = Traits<lb::D3Q15, lb::GuoForcingLBGK, lb::Normal,
                lb::DefaultStreamer, lb::DefaultWallStreamer, lb::DefaultInletStreamer, lb::DefaultOutletStreamer,
                stencil::FourPoint>;

      auto BuildTemplateContainer = [] (CellContainer const &cellContainer) ->  std::shared_ptr<TemplateCellContainer>
	{
	  auto templates = std::make_shared<TemplateCellContainer>();
	  for (auto& cell : cellContainer) {
	    templates->emplace(cell->GetTemplateName(), cell->clone());
	  }
	  return templates;
	};

      auto timers = std::make_unique<reporting::Timers>();
      
      SECTION("testCell2FluidWithoutCells") {
        CellContainer cells;
        CellArmy<Traits> army(*latDat, cells, BuildTemplateContainer(cells), *timers, cutoff);
        army.SetCell2Cell(/* intensity */1e0, /* cutoff */0.5);
        army.SetCell2Wall(/* intensity */1e0, /* cutoff */0.5);
        army.Cell2FluidInteractions();
      }

      SECTION("testCell2Fluid") {
	// Fixture all pairs far from one another
	auto cells = TwoPancakeSamosas<FakeCell>(cutoff);
	assert(cells.size() == 2);
	assert(std::dynamic_pointer_cast<FakeCell>( (*cells.begin()))->nbcalls == 0);
	assert(std::dynamic_pointer_cast<FakeCell>( (*std::next(cells.begin())))->nbcalls == 0);

	helpers::ZeroOutFOld(latDat.get());
	helpers::ZeroOutForces(latDat.get());

	CellArmy<Traits> army(*latDat, cells, BuildTemplateContainer(cells), *timers, cutoff);
	army.SetCell2Cell(/* intensity */1e0, /* cutoff */0.5);
	army.SetCell2Wall(/* intensity */1e0, /* cutoff */0.5);
	army.Cell2FluidInteractions();

	REQUIRE(std::dynamic_pointer_cast<FakeCell>( (*cells.begin()))->nbcalls == 1);
	REQUIRE(std::dynamic_pointer_cast<FakeCell>( (*std::next(cells.begin())))->nbcalls
		== 1);

	for (site_t i(0); i < latDat->GetDomain().GetLocalFluidSiteCount(); ++i) {
	  REQUIRE(latDat->GetSite(i).GetForce() == ApproxVector<LatticeForce>(0.0));
	}

	LatticePosition const n0(15 - 0.1, 15.5, 15.5);
	LatticePosition const n1(15 + 0.2, 15.5, 15.5);
	(*cells.begin())->GetVertices().front() = n0;
	(*std::next(cells.begin()))->GetVertices().front() = n1;
	army.updateDNC();
	army.Cell2FluidInteractions();

	REQUIRE(std::dynamic_pointer_cast<FakeCell>( (*cells.begin()))->nbcalls == 2);
	REQUIRE(std::dynamic_pointer_cast<FakeCell>( (*std::next(cells.begin())))->nbcalls
		== 2);
	REQUIRE(latDat->GetSite(15, 15, 15).GetForce() != ApproxVector<LatticeForce>(0.0));
      }

      SECTION("testCellInsertion") {
        auto cell = std::make_shared<FakeCell>(pancakeSamosa());
        // Shift cell to be contained in flow domain
        *cell += LatticePosition(2, 2, 2);

        helpers::ZeroOutFOld(latDat.get());
        helpers::ZeroOutForces(latDat.get());

        auto cells = CellContainer();
        CellArmy<Traits> army(*latDat, cells, BuildTemplateContainer(cells), *timers, cutoff);
        int called = 0;
        //using CellInserter = std::function<void(CellContainer::value_type)>;
        auto callback = [cell, &called](CellInserter const& inserter)
        {
          ++called;
          inserter(cell);
        };
        army.SetCellInsertion(callback);

        army.CallCellInsertion();
	REQUIRE(1 == called);
        REQUIRE(3ul == army.GetDNC().size());
        REQUIRE(1ul == army.GetCells().size());
      }

      SECTION("testFluid2Cell") {
        // Checks that positions of cells are updated. Does not check that attendant
        // DNC is.
        auto cells = TwoPancakeSamosas<FakeCell>(cutoff);
        auto const orig = TwoPancakeSamosas<FakeCell>(cutoff);
        auto const normal = Facet( (*cells.begin())->GetVertices(),
                                  (*cells.begin())->GetFacets()[0]).normal();

        LatticePosition gradient;
        Dimensionless non_neg_pop;
        std::function<Dimensionless(LatticeVelocity const &)> linear, linear_inv;
        std::tie(non_neg_pop, gradient, linear, linear_inv) = helpers::makeLinearProfile(cubeSizeWithHalo,
                                                                                         latDat.get(),
                                                                                         normal);

        CellArmy<Traits> army(*latDat, cells, BuildTemplateContainer(cells), *timers, cutoff);
        army.Fluid2CellInteractions();

        for (size_t i(0); i < cells.size(); ++i)
        {
          auto const disp = (*std::next(cells.begin(), i))->GetVertices().front()
              - (*std::next(orig.begin(), i))->GetVertices().front();
          auto i_nodeA = (*std::next(cells.begin(), i))->GetVertices().begin();
          auto i_nodeB = (*std::next(orig.begin(), i))->GetVertices().begin();
          auto const i_end = (*std::next(cells.begin(), i))->GetVertices().end();

          for (; i_nodeA != i_end; ++i_nodeA, ++i_nodeB)
          {
            REQUIRE((*i_nodeA - *i_nodeB) == ApproxV(disp));
          }
        }
      }

      SECTION("testCellOutput") {
        auto cell = std::make_shared<FakeCell>(tetrahedron());
        // Shift cell to be contained in flow domain
        *cell += LatticePosition(2, 2, 2);

        MeshData::Vertices::value_type barycentre;
        CellChangeListener callback =
            [&barycentre](const CellContainer & container)
            {
              barycentre = (*(container.begin()))->GetBarycentre();
            };

        CellContainer intel;
        intel.insert(cell);
        CellArmy<Traits> army(*latDat, intel, BuildTemplateContainer(intel), *timers, cutoff);
        army.AddCellChangeListener(callback);

        army.NotifyCellChangeListeners();
	REQUIRE(barycentre == cell->GetBarycentre());
      }

      SECTION("testCellRemoval") {
        // The flow extension takes the z > CubeSize()/2 half of the cube
	int len = cubeSizeWithHalo / 2;
        FlowExtension const outlet(util::Vector3D<Dimensionless>(0, 0, 1),
                                   LatticePosition(len, len, len),
                                   len,
                                   len,
                                   cubeSizeWithHalo / 4);
        auto cell = std::make_shared<FakeCell>(pancakeSamosa());
        // Shift cell to be contained in flow domain
        *cell += LatticePosition(2, 2, 2);

        helpers::ZeroOutFOld(latDat.get());
        helpers::ZeroOutForces(latDat.get());

        CellContainer intel;
        intel.insert(cell);
        CellArmy<Traits> army(*latDat, intel, BuildTemplateContainer(intel), *timers, cutoff);
        army.SetOutlets(std::vector<FlowExtension>(1, outlet));

        // Check status before attempting to remove cell that should *not* be removed
        REQUIRE(3ul == army.GetDNC().size());
        REQUIRE(1ul == army.GetCells().size());

        // Check status after attempting to remove cell that should *not* be removed
        army.CellRemoval();
        REQUIRE(3ul == army.GetDNC().size());
        REQUIRE(1ul == army.GetCells().size());

        // Check status after attempting to remove cell that should be removed.
        // Move the cell beyond the fading section of the flow extension.
        *cell += LatticePosition(outlet.origin + outlet.normal * outlet.fadeLength * 1.1);
        army.CellRemoval();
        REQUIRE(0ul == army.GetDNC().size());
        REQUIRE(0ul == army.GetCells().size());
      }

    }

}

