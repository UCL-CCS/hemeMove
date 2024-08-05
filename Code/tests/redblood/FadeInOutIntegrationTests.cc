// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <memory>
#include <catch2/catch.hpp>

#include "SimulationController.h"
#include "configuration/SimBuilder.h"
#include "lb/lattices/D3Q19.h"
#include "Traits.h"
#include "redblood/Mesh.h"
#include "redblood/Cell.h"
#include "redblood/CellController.h"
#include "tests/helpers/FolderTestFixture.h"

namespace hemelb::tests
{
    using namespace redblood;
    TEST_CASE_METHOD(helpers::FolderTestFixture,
                     "FadeInOutIntegrationTests",
                     "[redblood][.long]") {
        using Traits = Traits<lb::D3Q19, lb::GuoForcingLBGK>;
        using CellControl = CellController<Traits>;

      CopyResourceToTempdir("large_cylinder_rbc.xml");
      CopyResourceToTempdir("large_cylinder.gmy");
      CopyResourceToTempdir("red_blood_cell.txt");

      ModifyXMLInput("large_cylinder_rbc.xml", { "simulation", "steps", "value" }, 22000);
      ModifyXMLInput("large_cylinder_rbc.xml",
		     { "redbloodcells", "controller", "stencil" },
		     "two");
      ModifyXMLInput("large_cylinder_rbc.xml", { "outlets",
	    "outlet",
	    "flowextension",
	    "length",
	    "value" },
	40e-6);

      constexpr int argc = 7;
      char const * argv[argc] = {
	"hemelb", "-in", "large_cylinder_rbc.xml",
      };

      auto options = configuration::CommandLine{argc, argv};
      auto sim = configuration::SimBuilder::CreateSim<Traits>(options, Comms());
     
      SECTION("testIntegration") {
	auto const & converter = sim->GetUnitConverter();
	auto const volumeFactor = std::pow(converter.ConvertToLatticeUnits("m", 1e0), -3)
	  * 1e12;
	bool didDropCell = false;
	auto checkDidDropCell = [&didDropCell](const CellContainer & cells)
	  {
	    didDropCell |= not cells.empty();
	  };
	auto checkVolume = [volumeFactor](const hemelb::redblood::CellContainer & cells) {
	  static LatticeVolume expected = -1e0;
	  if(cells.empty()) {
	    return;
	  }
	  REQUIRE(1 == int(cells.size()));
	  auto cell = *cells.begin();
	  auto const volume = cell->GetVolume() * volumeFactor;
	  if (expected < 0e0) {
	    expected = volume;
	  }
	  REQUIRE(Approx(expected).margin(1e-7) == volume);
	};
	auto checkPosition = []( const hemelb::redblood::CellContainer & cells) {
	  if(cells.empty()) {
	    return;
	  }
	  REQUIRE(1 == int(cells.size()));
	  auto cell = *cells.begin();
	  static LatticePosition first, current, tenth;
	  static int iter = 0;
	  if(iter == 0) {
	    first = cell->GetBarycentre();
	    current = first;
	    tenth = first;
	  }
	  auto const position = cell->GetBarycentre();
	  REQUIRE(Approx(first.x()).margin(1e-8) == position.x());
	  REQUIRE(Approx(first.y()).margin(1e-8) == position.y());
	  REQUIRE(current.z() <= position.z());
	  current = position;
	  ++iter;
	  if(iter % 10 == 0) {
	    REQUIRE(tenth.z() < position.z());
	    tenth = position;
	  }
	};
	int iter = 0;
	auto iterate = [&iter](const hemelb::redblood::CellContainer&) {
	  ++iter;
	};

	REQUIRE(sim);
	auto controller = std::static_pointer_cast<CellControl>(sim->GetCellController());
	REQUIRE(controller);
	controller->AddCellChangeListener(checkVolume);
	controller->AddCellChangeListener(checkPosition);
	controller->AddCellChangeListener(iterate);
	controller->AddCellChangeListener(checkDidDropCell);

	// keep those lambdas inline to avoid unused function warning when commented out.
	controller->AddCellChangeListener([]( const hemelb::redblood::CellContainer & cells) {
	    static int iter = -1;
	    ++iter;
	    if(cells.empty()) {
	      return;
	    }
	    auto cell = *cells.begin();
	    auto const tag = cell->GetTag();
	    auto const b = cell->GetBarycentre();
	    auto const v = cell->GetVolume();
	    auto const e = (*cell)();
	    HEMELB_CAPTURE5(iter, tag, b, v, e);
	  });

	// controller->AddCellChangeListener(
	//     [&converter](const hemelb::redblood::CellContainer &cells)
	//     {
	//       static int iter = 0;
	//       if(cells.empty())
	//       {
	//         return;
	//       }
	//       auto cell = *cells.begin();
	//       if(iter% 1000 == 0)
	//       {
	//         std::ostringstream sstr;
	//         sstr << "/tmp/cell-" << cell->GetTag() << "_" << iter<< ".vtp";
	//         writeVTKMesh(sstr.str(), cell, converter);
	//       }
	//       ++iter;
	//     }
	// );

	// run the simulation
	sim->RunSimulation();
	sim->Finalise();

	AssertPresent("results/report.txt");
	AssertPresent("results/report.xml");
	REQUIRE(iter > 0);
	REQUIRE(didDropCell);
	REQUIRE(controller->GetCells().empty());
      }
    }

}

