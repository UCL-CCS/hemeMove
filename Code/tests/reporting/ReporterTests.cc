// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <iomanip>

#include <catch2/catch.hpp>
#include <ctemplate/template.h>

#include "build_info.h"
#include "lb/IncompressibilityChecker.h"
#include "lb/IncompressibilityChecker.hpp"
#include "reporting/BuildInfo.h"
#include "reporting/Reporter.h"
#include "reporting/Timers.h"
#include "reporting/Timers.hpp"
#include "util/Iterator.h"

#include "tests/helpers/FourCubeLatticeData.h"
#include "tests/helpers/HasCommsTestFixture.h"
#include "tests/lb/BroadcastMocks.h"
#include "tests/lb/LbTestsHelper.h"
#include "tests/reporting/Mocks.h"

namespace hemelb::tests
{
    using namespace hemelb::reporting;

    using TimersMock = TimersBase<ClockMock>;
    using IncompressibilityCheckerMock = lb::IncompressibilityChecker<net::BroadcastMockRootNode>;

    TEST_CASE_METHOD(helpers::HasCommsTestFixture, "ReporterTests") {

        // We need two sets of timers because the incompressibility
        // checker is not templated over timing policy.
        auto mockTimers = TimersMock();
        auto mockComms = MPICommsMock();
        auto realTimers = reporting::Timers();
        auto buildInfo = reporting::BuildInfo();
        auto state = lb::SimulationState(0.0001, 1000);
        auto net = net::Net(Comms());
        auto latticeData = std::unique_ptr<FourCubeLatticeData>(FourCubeLatticeData::Create(Comms(), 6, 5)); // The 5 here is to match the topology size in the MPICommsMock
        auto& domain = latticeData->GetDomain();
        LbTestsHelper::InitialiseAnisotropicTestData<lb::D3Q15>(*latticeData);
        latticeData->SwapOldAndNew(); //Needed since InitialiseAnisotropicTestData only initialises FOld
        auto cache = lb::MacroscopicPropertyCache(state, domain);
        cache.densityCache.SetRefreshFlag();
        LbTestsHelper::UpdatePropertyCache<lb::D3Q15>(*latticeData, cache, state);
        auto incompChecker = IncompressibilityCheckerMock(&domain, &net, &state, cache, realTimers, 10.0);
        auto reporter = Reporter("mock_path", "exampleinputfile");
        reporter.AddReportable(&incompChecker);
        reporter.AddReportable(&mockTimers);
        reporter.AddReportable(&state);
        reporter.AddReportable(&domain);
        reporter.AddReportable(&buildInfo);

        auto AssertTemplate = [&](const std::string &expectation, const std::string &ttemplate) {
            ctemplate::StringToTemplateCache("TestFor" + ttemplate, ttemplate, ctemplate::DO_NOT_STRIP);
            std::string result;
            REQUIRE(ctemplate::ExpandTemplate("TestFor" + ttemplate,
                                              ctemplate::DO_NOT_STRIP,
                                              reporter.GetDictionary().GetRaw(),
                                              &result));
            REQUIRE(expectation == result);
        };

        auto AssertValue = [&](const std::string & expectation, const std::string &symbol) {
            AssertTemplate(expectation, "{{" + symbol + "}}");
        };

        auto CheckTimingsTable = [&]() {
            std::stringstream expectation;
            expectation << std::setprecision(3);
            for (auto [row, t]: util::enumerate(mockTimers)) {
                expectation << "N" << t.Description() << "L" << row * 10.0 << "MI" << row * 15.0 << "ME"
                            << row * 2.0 << "MA" << row * 5.0 << " " << std::flush;
            }
            AssertTemplate(expectation.str(), "{{#TIMER}}N{{NAME}}L{{LOCAL}}MI{{MIN}}ME{{MEAN}}MA{{MAX}} {{/TIMER}}");
        };

        SECTION("TestInit") {
            AssertValue("exampleinputfile", "CONFIG");
        }



        SECTION("TestMainReport") {
            // Mock up some timings
            for (auto [i, t]: util::enumerate(mockTimers)) {
                for (unsigned int j = 0; j < i; j++) {
                    t.Start();
                    t.Stop();
                }
            }
            mockTimers.Reduce(mockComms); // invoke the Timers MPI mock
            for (unsigned int step = 0; step < 1000; step++) {
                state.Increment();
            }
            REQUIRE(1000lu == state.GetTimeStep());
            reporter.FillDictionary();

            CheckTimingsTable();
            AssertTemplate("", "{{#UNSTABLE}} unstable{{/UNSTABLE}}");
            AssertTemplate("R0S64 R1S1000 R2S2000 R3S3000 R4S4000 ",
                           "{{#PROCESSOR}}R{{RANK}}S{{SITES}} {{/PROCESSOR}}");
            AssertTemplate(build_info::REVISION_HASH.str(), "{{#BUILD}}{{REVISION}}{{/BUILD}}");
            AssertTemplate(build_info::BUILD_TIME.str(), "{{#BUILD}}{{TIME}}{{/BUILD}}");
            AssertValue("0.000100", "TIME_STEP_LENGTH");
            AssertValue("1000", "TOTAL_TIME_STEPS");
            AssertValue("1000", "STEPS");
            AssertValue("64", "SITES");
            AssertValue("1", "BLOCKS");
            AssertValue("216", "SITESPERBLOCK");
        }

    }
}
