// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <catch2/catch.hpp>

#include "resources/Resource.h"
#include "multiscale/mpwide/MPWideIntercommunicator.h"
#include "multiscale/MultiscaleSimulationController.h"

#include "tests/helpers/FolderTestFixture.h"
#include "tests/multiscale/MockIntercommunicand.h"
#include "tests/multiscale/mpwide/MockMPWide.h"
#include "tests/multiscale/mpwide/IntercommunicatingHemeLB.h"

namespace hemelb
{
  namespace tests
  {
    using namespace multiscale;
    
    TEST_CASE_METHOD(helpers::FolderTestFixture, "MPWideIntercommunicatorTests") {
	InterCommunicatingHemeLB<MPWideIntercommunicator> *mockheme;
	std::map<std::string, double> *pbuffer;
	std::map<std::string, bool> *LBorchestration;

	
	pbuffer = new std::map<std::string, double>();

	LBorchestration = new std::map<std::string, bool>();
	(*LBorchestration)["boundary1_pressure"] = false;
	(*LBorchestration)["boundary2_pressure"] = false;
	(*LBorchestration)["boundary1_velocity"] = true;
	(*LBorchestration)["boundary2_velocity"] = true;

	CopyResourceToTempdir("MPWSettings.cfg");
	std::string configPath = "MPWSettings.cfg";
	mockheme = new InterCommunicatingHemeLB<MPWideIntercommunicator>(25.0,
									 0.2,
									 *pbuffer,
									 *LBorchestration,
									 configPath);
	SECTION("testMPWidePresent") {
	  std::string host = "localhost";
	  std::cout << "IP address for localhost is: " << MPW_DNSResolve(const_cast<char*>(host.c_str()))
		    << std::endl;
	  std::cout << "MPWide is present." << std::endl;
	}

	SECTION("testMPWideInit") {
	  // TODO This test needs writing.
	}

	SECTION("testMPWideApplication") {
	  int constexpr argc = 3;
	  const char* argv[argc] = {
	    "hemelb",
	    "-in",
	    "four_cube_multiscale.xml",
	  };

	  CopyResourceToTempdir("four_cube_multiscale.xml");
	  CopyResourceToTempdir("four_cube.gmy");
	  configuration::CommandLine options(argc, argv);
	  std::string configPath = "../../../config_files/MPWSettings.cfg";
	  MPWideIntercommunicator intercomms(Comms().OnIORank(), *pbuffer, *LBorchestration, configPath);

	  MultiscaleSimulationController<MPWideIntercommunicator> heme(options, Comms(), intercomms);
	  // Mock out the behaviour of the simulation iteration, but with the other model linked in.
	  //std::cout << "HemeLB about to be run..." << std::endl;
	  while (heme.GetState()->GetTime() < 20.0) {
	    heme.DoTimeStep();
	    //std::cout << "Step taken, going to incrementSharedTime." << std::endl;
	    intercomms.UnitTestIncrementSharedTime(); //simple hack func that mocks a 1.0 increase in the 'other' simulation.
	  }
	  heme.Finalise();
	  REQUIRE(heme.GetState()->GetTime() == Approx(20.0));
	  //CPPUNIT_ASSERT_DOUBLES_EQUAL(zerod->currentTime, 20.5, 1e-6); // does one more step, where it sets the shared time.
	}

	delete mockheme;
	delete pbuffer;
	delete LBorchestration;
      }

  } //tests
} //hemelb
