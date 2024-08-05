// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#include <memory>

#include <catch2/catch.hpp>

#include "configuration/SimBuilder.h"

#include "tests/helpers/FolderTestFixture.h"
#include "tests/helpers/LaddFail.h"

namespace fs = std::filesystem;
namespace hemelb::tests
{

    TEST_CASE_METHOD(helpers::FolderTestFixture, "SimulationController") {
      const int argc = 3;
      const char* argv[] = {
	"hemelb",
	"-in",
	"four_cube.xml",
      };

      CopyResourceToTempdir("four_cube.xml");
      CopyResourceToTempdir("four_cube.gmy");

      auto options = std::make_unique<hemelb::configuration::CommandLine>(argc, argv);
      auto controller = configuration::SimBuilder::CreateSim<Traits<>>(*options, Comms());

      SECTION("Running a simulation creates outputs") {
	// TODO: This test is fatal if run with LADDIOLET. See ticket #605.
	LADD_FAIL();
	controller->RunSimulation();
    for (auto&& p: {"results/report.txt",
                    "results/report.xml",
                    "results/Extracted/wholegeometryvelocityandstress.dat",
                    "results/Extracted/centrelinepressure.dat",
                    "results/Extracted/centrelineshearrate.dat",
                    "results/Extracted/surfaceshearstress.dat",
                    "results/Extracted/surfacetraction.dat"}) {
        INFO(p);
        REQUIRE(fs::exists(p));
    }
      }
    }
}
