// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <catch2/catch.hpp>

#include "reporting/Timers.h"
#include "reporting/Timers.hpp"

#include "util/Iterator.h"
#include "tests/reporting/Mocks.h"
#include "tests/helpers/HasCommsTestFixture.h"

namespace hemelb::tests
{
    using namespace hemelb::reporting;

    TEST_CASE("TimerTests") {
        auto timer = TimerBase<ClockMock>();

        SECTION("TestInitialization") {
            // should default to zero time on startup
            REQUIRE(0.0 == timer.Get());
        }

        SECTION("TestStartStop") {
            timer.Start(); // clock mock at 10.0
            timer.Stop(); // clock mock at 20.0
            REQUIRE(Approx(10.0) == timer.Get());
        }

        SECTION("TestMultipleStartStop") {
            timer.Start(); // clock mock at 10.0
            timer.Stop(); // clock mock at 20.0
            timer.Start(); // clock mock at 30.0
            timer.Stop(); // clock mock at 40.0
            REQUIRE(Approx(20.0) == timer.Get());
        }

        SECTION("TestSetTime") {
            timer.Set(15.0);
            REQUIRE(Approx(15.0) == timer.Get());
            timer.Start(); // clock mock at 10.0
            timer.Stop(); // clock mock at 20.0
            REQUIRE(Approx(25.0) == timer.Get());
        }
    }

    TEST_CASE("TimersTests") {
        auto timers = TimersBase<ClockMock>();

        SECTION("TestInitialization") {
            REQUIRE(Approx(0.0) == timers.total().Get());
            for (auto& t: timers) {
                REQUIRE(Approx(0.0) == t.Get());
            }
        }

        SECTION("TestTimersSeparate") {
            for (auto [i, t]: util::enumerate(timers)) {
                for (unsigned int j = 0; j < i; j++) {
                    t.Start();
                    t.Stop();
                }
            }
            for (auto [i, t]: util::enumerate(timers)) {
                REQUIRE(Approx(10.0 * i) == t.Get());
            }
        }

        SECTION("TestReduce") {
            for (auto [i, t]: util::enumerate(timers)) {
                for (unsigned int j = 0; j < i; j++) {
                    t.Start();
                    t.Stop();
                }
            }

            timers.Reduce(MPICommsMock()); // trigger the mock
            for (unsigned int i = 0; i < Timers::numberOfTimers; i++) {
                REQUIRE(Approx(i * 5.0) == timers.Maxes()[i]);
                REQUIRE(Approx(i * 2.0) == timers.Means()[i]);
                REQUIRE(Approx(i * 15.0) == timers.Mins()[i]);
            }
        }

    }
}
