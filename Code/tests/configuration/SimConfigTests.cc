// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <memory>

#include <catch2/catch.hpp>

#include "configuration/SimConfig.h"
#include "resources/Resource.h"
#include "tests/helpers/FolderTestFixture.h"
#include "tests/helpers/LaddFail.h"

namespace hemelb::tests
{
    using namespace configuration;
    namespace
    {
      // For section XMLFileContent
      struct CfgChecker {
	using result_type = bool;
	template<typename T>
	bool operator()(T) const {
	  return false;
	}
	bool operator()(const EquilibriumIC& eqIC) const {
	  return 80.0 == eqIC.p_Pa;
	}
      };
    }
    TEST_CASE_METHOD(helpers::FolderTestFixture, "SimConfig") {

      SECTION("0_2_0_Read") {
	LADD_FAIL();
	// smoke test the configuration as having loaded OK
	auto config = SimConfig::New(resources::Resource("config0_2_0.xml").Path());
	REQUIRE(3000lu == config.GetTotalTimeSteps());
	REQUIRE(0.0001 == config.GetTimeStepLength());
    auto inlet = std::get_if<configuration::CosinePressureIoletConfig>(&config.GetInlets()[0]);
	//auto inlet = util::clone_dynamic_cast<lb::iolets::InOutLetCosine>(config->GetInlets()[0]);

	REQUIRE(inlet != nullptr);
	REQUIRE(Approx(0.6) == inlet->period_s);

	// Check that in the absence of the <monitoring> XML element things get initiliased properly
	auto& monConfig = config.GetMonitoringConfiguration();
	REQUIRE(!monConfig.doConvergenceCheck);
	REQUIRE(!monConfig.doIncompressibilityCheck);
	REQUIRE(!monConfig.convergenceTerminate);
	REQUIRE(Approx(0.0) == monConfig.convergenceRelativeTolerance);
      }

      SECTION("0_2_1_Read") {
	LADD_FAIL();
	// smoke test the configuration as having loaded OK
	auto config = SimConfig::New(resources::Resource("config.xml").Path());
	REQUIRE(3000lu == config.GetTotalTimeSteps());
	REQUIRE(0.0001 == config.GetTimeStepLength());
	auto inlet = std::get_if<configuration::CosinePressureIoletConfig>(&config.GetInlets()[0]);
	REQUIRE(inlet != nullptr);
	REQUIRE(Approx(0.6) == inlet->period_s);

	auto& monConfig = config.GetMonitoringConfiguration();
	REQUIRE(monConfig.doConvergenceCheck);
	REQUIRE(monConfig.doIncompressibilityCheck);
	REQUIRE(monConfig.convergenceTerminate);
	REQUIRE(1e-9 == monConfig.convergenceRelativeTolerance);
	REQUIRE(std::holds_alternative<extraction::source::Velocity>(monConfig.convergenceVariable));
	REQUIRE(1 == monConfig.convergenceReferenceValue); // 1 m/s
      }
      
      SECTION("XMLFileContent") {
	LADD_FAIL();
	//Round trip the config twice.
	CopyResourceToTempdir("config.xml");
	auto config = SimConfig::New("config.xml");

	auto& ICconfig = config.GetInitialCondition();
	REQUIRE(std::visit(CfgChecker{}, ICconfig));
      }

    }
}
