// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <catch2/catch.hpp>
#include <tinyxml2.h>

#include "Traits.h"
#include "SimulationController.h"
#include "configuration/SimBuilder.h"
#include "redblood/Cell.h"
#include "redblood/CellController.h"
#include "tests/helpers/LatticeDataAccess.h"
#include "tests/helpers/FolderTestFixture.h"

namespace hemelb::tests
{
    using namespace redblood;
    class CellIntegrationTests : public helpers::FolderTestFixture
    {
      using MyTraits = Traits<lb::DefaultLattice, lb::GuoForcingLBGK>;
      using CellControll = CellController<MyTraits>;

    public:
      CellIntegrationTests() : FolderTestFixture(), timings() {
	CopyResourceToTempdir("red_blood_cell.txt");
	tinyxml2::XMLDocument doc;
    doc.LoadFile(resources::Resource("large_cylinder.xml").Path().c_str());
	CopyResourceToTempdir("large_cylinder.xml");
    ModifyXMLInput("large_cylinder.xml", {"simulation", "steps", "value"}, 8);
	CopyResourceToTempdir("large_cylinder.gmy");
	options = std::make_unique<configuration::CommandLine>(argc, argv);

	auto cell = std::make_shared<Cell>(icoSphere(4));
	templates = std::make_unique<TemplateCellContainer>();
	(*templates)["icosphere"] = cell;
	cell->moduli = Cell::Moduli(1e-6, 1e-6, 1e-6, 1e-6);
	cells.insert(cell);

	//timings = std::make_unique<reporting::Timers>(Comms());
	simController = configuration::SimBuilder::CreateSim<MyTraits>(*options, Comms());
	helpers::LatticeDataAccess(&simController->GetFieldData()).ZeroOutForces();
      }

      ~CellIntegrationTests() {
	simController->Finalise();
      }

      // No errors when interpolation/spreading hits nodes out of bounds
      void testCellOutOfBounds()
      {
          auto& fd = simController->GetFieldData();
          auto& dom = fd.GetDomain();
	(*cells.begin())->operator+=(dom.GetGlobalSiteMins() * 2.0);
	auto controller = std::make_shared<CellControll>(
							 fd,
							 cells,
							 templates,
							 timings
							 );
	simController->RegisterActor(*controller, 1);
	simController->RunSimulation();
	AssertPresent("results/report.txt");
	AssertPresent("results/report.xml");
      }

      // Check that the particles move and result in some force acting on the fluid
      void testIntegration()
      {
	// setup cell position
	auto& fieldData = simController->GetFieldData();
    auto& dom = fieldData.GetDomain();
	auto const mid = LatticePosition(dom.GetGlobalSiteMaxes()
                                     + dom.GetGlobalSiteMins()) * 0.5;
	(**cells.begin()) += mid - (*cells.begin())->GetBarycentre();
	(**cells.begin()) += LatticePosition(0, 0, 8 - mid.z());
	(**cells.begin()) *= 5.0;
	auto controller = std::make_shared<CellControll>(fieldData, cells, templates, timings);
	auto const barycentre = (*cells.begin())->GetBarycentre();

	// run
	simController->RegisterActor(*controller, 1);
	simController->RunSimulation();

	// check position of cell has changed
	auto const moved = (*cells.begin())->GetBarycentre();
	REQUIRE(Approx(barycentre.x()).margin(1e-12) == moved.x());
	REQUIRE(Approx(barycentre.y()).margin(1e-12) == moved.y());
	REQUIRE(std::abs(barycentre.z() - moved.z()) > 1e-8);

	// Check there is force on one of the lattice site near a
	// node node position is guessed at from geometry. This
	// truncates.
	auto const nodepos = LatticeVector{mid + LatticePosition(0, 0, 8 - 5 - mid.z())};
	auto const force = fieldData.GetSite(nodepos).GetForce();
	REQUIRE(std::abs(force.z()) > 1e-4);

	AssertPresent("results/report.txt");
	AssertPresent("results/report.xml");
      }

      // Check that the particles move and result in some force acting on the fluid
      void testIntegrationWithoutCells()
      {
	// setup cell position
	CellContainer empty;
	auto empty_tmpl = std::make_shared<TemplateCellContainer>();
	auto controller = std::make_shared<CellControll>(
            simController->GetFieldData(),
            empty, empty_tmpl, timings);

	// run
	simController->RegisterActor(*controller, 1);
	simController->RunSimulation();

	AssertPresent("results/report.txt");
	AssertPresent("results/report.xml");
      }

    private:
      std::unique_ptr<SimulationController> simController;
      std::unique_ptr<configuration::CommandLine> options;
      CellContainer cells;
      std::shared_ptr<TemplateCellContainer> templates;
      reporting::Timers timings;
      static constexpr int argc = 3;
      static char const* const argv[argc];
    };
    char const* const CellIntegrationTests::argv[]  = {"hemelb", "-in", "large_cylinder.xml"};

    METHOD_AS_TEST_CASE(CellIntegrationTests::testCellOutOfBounds,
			"No errors when interpolation/spreading hits nodes out of bounds",
			"[redblood][.long]");
    METHOD_AS_TEST_CASE(CellIntegrationTests::testIntegration,
			"Check that the particles move and result in some force acting on the fluid",
			"[redblood][.long]");
    METHOD_AS_TEST_CASE(CellIntegrationTests::testIntegrationWithoutCells,
			"Test things work without any cells",
			"[redblood][.long]");

}
