// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <iostream>

#include <catch2/catch.hpp>

#include "resources/Resource.h"
#include "multiscale/MultiscaleSimulationController.h"

#include "tests/helpers/LaddFail.h"
#include "tests/helpers/FolderTestFixture.h"
#include "tests/multiscale/MockIntercommunicand.h"
#include "tests/multiscale/MockIntercommunicator.h"

namespace hemelb::tests
{
    using namespace multiscale;
    /***
     * Mock intercommunicating entity which looks a bit like a HemeLB conceptually
     * It has an input, and an output, and the flow rate depends on the difference in pressures.
     */
    template<class IntercommuniatorImplementation>
    class MockHemeLB
    {
    public:
        MockHemeLB(double spaceResolution, double timeResolution,
                   std::map<std::string, double> & buffer,
                   std::map<std::string, bool> &orchestration) :
                inlet(1.0, 0.1), outlet(-1.0, 0.1), inOutLetType("inoutlet"),
                intercomms(buffer, orchestration), timeResolution(timeResolution),
                spaceResolution(spaceResolution), currentTime(0)
        {
            // The intercommunicators have a shared buffer which represents imaginary communication
            inOutLetType.template RegisterSharedValue<double>("pressure");
            inOutLetType.template RegisterSharedValue<double>("velocity");
            intercomms.RegisterIntercommunicand(inOutLetType, inlet, "boundary1");
            intercomms.RegisterIntercommunicand(inOutLetType, outlet, "boundary2");
            intercomms.ShareInitialConditions();
        }
        MockIntercommunicand inlet;
        MockIntercommunicand outlet;
        typename IntercommuniatorImplementation::IntercommunicandTypeT inOutLetType;
        IntercommuniatorImplementation intercomms;
        double timeResolution;
        double spaceResolution;
        double currentTime;
        void DoLB()
        {
            double resistance = 10.0;
            double velocity = (inlet.GetPressure() - outlet.GetPressure()) / resistance;
            inlet.SetVelocity(velocity);
            outlet.SetVelocity(velocity);
        }

        void Simulate()
        {
            if (intercomms.DoMultiscale(currentTime))
            {
                log::Logger::Log<log::Debug, log::OnePerCore>("DoLB() MH currentTime: %f time Resolution: %f P in/out: %f %f",
                                                                                      currentTime,
                                                                                      timeResolution,
                                                                                      inlet.GetPressure(),
                                                                                      outlet.GetPressure());
                DoLB();
                currentTime += timeResolution;
            }
        }
    };

    /***
     * Mock intercommunicating entity which is a tank of water.
     * It has a pressure difference which drops depending on the flow rate.
     */
    template<class IntercommuniatorImplementation>
    class Mock0DModel
    {
    public:
        Mock0DModel(double spaceResolution, double timeResolution,
                    std::map<std::string, double> & buffer,
                    std::map<std::string, bool> &orchestration) :
                inlet(-1.0, 0.1), outlet(1.0, 0.1), inOutLetType("inoutlet"),
                intercomms(buffer, orchestration), timeResolution(timeResolution),
                spaceResolution(spaceResolution), currentTime(0)
        {
            // The intercommunicators have a shared buffer which represents imaginary communication
            inOutLetType.template RegisterSharedValue<double>("pressure");
            inOutLetType.template RegisterSharedValue<double>("velocity");
            intercomms.RegisterIntercommunicand(inOutLetType, inlet, "boundary2");
            intercomms.RegisterIntercommunicand(inOutLetType, outlet, "boundary1");
            intercomms.ShareInitialConditions();
        }
        MockIntercommunicand inlet;
        MockIntercommunicand outlet;
        typename IntercommuniatorImplementation::IntercommunicandTypeT inOutLetType;
        IntercommuniatorImplementation intercomms;
        double timeResolution;
        double spaceResolution;
        double currentTime;
        double GetOutletPressure()
        {
            return outlet.GetPressure();
        }
        double GetTime()
        {
            return currentTime;
        }
        void Do1D()
        {
            double capacitance = 10.0;
            double deltap = -1.0 * inlet.GetVelocity() / capacitance;
            log::Logger::Log<log::Debug, log::OnePerCore>("0D (deltap, inlet vel): %f, %f",
                                                          deltap,
                                                          inlet.GetVelocity());
            outlet.SetPressure(outlet.GetPressure() + deltap);
        }

        void Simulate()
        {
            if (intercomms.DoMultiscale(currentTime))
            {
                Do1D();
                currentTime += timeResolution;
                log::Logger::Log<log::Debug, log::OnePerCore>("0D: %f, %f",
                                                              currentTime,
                                                              GetOutletPressure());
            }
        }

    };

    // Useful for debugging to have this.
    std::ostream & operator <<(std::ostream & stream, std::map<std::string, double> buffer)
    {
        stream << " { ";
        for (std::map<std::string, double>::iterator entry = buffer.begin(); entry != buffer.end();
             entry++)
        {
            stream << entry->first << " : " << entry->second << " , ";
        }
        stream << " } ";
        stream << std::flush;
        return stream;
    }

    TEST_CASE_METHOD(public helpers::FolderTestFixture, "MockIntercommunicatorTests") {
        MockHemeLB<MockIntercommunicator> *mockheme;
        MultiscaleSimulationController<MockIntercommunicator> *heme;
        Mock0DModel<MockIntercommunicator> *zerod;
        std::map<std::string, double> *pbuffer;
        std::map<std::string, bool> *orchestrationLB;
        std::map<std::string, bool> *orchestrationOD;

        pbuffer = new std::map<std::string, double>();
        std::map<std::string, double> &buffer = *pbuffer;

        orchestrationOD = new std::map<std::string, bool>();
        orchestrationLB = new std::map<std::string, bool>();
        std::map<std::string, bool> &rorchestrationLB = *orchestrationLB;
        std::map<std::string, bool> &rorchestration0D = *orchestrationOD;
        rorchestrationLB["boundary1_pressure"] = false;
        rorchestrationLB["boundary2_pressure"] = false;
        rorchestrationLB["boundary1_velocity"] = true;
        rorchestrationLB["boundary2_velocity"] = true;
        rorchestration0D["boundary1_pressure"] = true;
        rorchestration0D["boundary2_pressure"] = true;
        rorchestration0D["boundary1_velocity"] = false;
        rorchestration0D["boundary2_velocity"] = false;

        mockheme = new MockHemeLB<MockIntercommunicator>(25.0, 0.2, buffer, rorchestrationLB);
        zerod = new Mock0DModel<MockIntercommunicator>(10.0, 0.5, buffer, rorchestration0D);

        SECTION("TestCRRun") {
            // In this example, we imagine we are coupling a mock hemelb with a zero-d capacitative model.
            // Our example intercommunicator has been hard coded with how to share the values.
            // In principle, though, this would be defined by whatever configuration file format the implementing intercommunicator wants to use.
            for (int steps = 0; steps < 100; steps++)
            {
                mockheme->Simulate();
                zerod->Simulate();
            }
            // The asserted value here is a result of discretised exponential decay.
            // I have not thought it worth determining this by a formula.
            REQUIRE(Approx(0.3258457641491) == zerod->GetOutletPressure());
            REQUIRE(Approx(mockheme->currentTime) == 20.0);
            REQUIRE(Approx(zerod->currentTime) == 20.5);
        }

        SECTION("TestCHemeRun") {
            //CPPUNIT_ASSERT_MESSAGE("This test is broken - see ticket #663", 1 == 0);

            // TODO: This test is fatal if run with LADDIOLET. See ticket #605.
            LADD_FAIL();
            int argc;
            const char* argv[7];
            argc = 3;
            argv[0] = "hemelb";
            argv[2] = "four_cube_multiscale.xml";
            argv[1] = "-in";
            CopyResourceToTempdir("four_cube_multiscale.xml");
            CopyResourceToTempdir("four_cube.gmy");
            configuration::CommandLine options(argc, argv);
            MockIntercommunicator intercomms(*pbuffer, *orchestrationLB);
            heme = new MultiscaleSimulationController<MockIntercommunicator>(options,
                                                                         Comms(),
                                                                         intercomms);
            // Mock out the behaviour of the simulation iteration, but with the other model linked in.
            while (heme->GetState().GetTime() < 20.0 || zerod->GetTime() < 20.0)
            {
                heme->DoTimeStep();
                zerod->Simulate();
            }
            // In advancing 100 time steps, at 0.2 s per time step, with a 0d model at 0.5s per time step
            // the 0d model with execute 100*2/5=40 times, plus one more step, the one where it communicates the previous step.
            // Each time, the pressure difference will drop by 0.1*0.1=0.01 Pa.
            // So the final pressure will be 81.0-41*0.01=80.59 Pa
            REQUIRE(Approx(0.59) == zerod->GetOutletPressure());
            heme->Finalise();
            REQUIRE(heme->GetState().GetTime() == Approx(20.0));
            REQUIRE(zerod->currentTime == Approx(20.5)); // does one more step, where it sets the shared time.
            delete heme;
        }

        delete zerod;
        delete mockheme;
        delete pbuffer;
        //delete orchestrationLB;
    }

}
