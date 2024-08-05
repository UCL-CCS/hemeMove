// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <memory>
#include <boost/uuid/uuid_io.hpp>
#include <catch2/catch.hpp>

#include "configuration/SimBuilder.h"
#include "Traits.h"
#include "redblood/Mesh.h"
#include "redblood/MeshIO.h"
#include "redblood/CellController.h"
#include "tests/helpers/FolderTestFixture.h"

#ifdef HEMELB_CALLGRIND
#include <valgrind/callgrind.h>
#endif

namespace hemelb::tests
{
    using namespace redblood;

    class MulticellBenchmarkTests : public helpers::FolderTestFixture
    {
        using MyTraits = Traits<lb::DefaultLattice, lb::GuoForcingLBGK>;
        using CellControl = hemelb::redblood::CellController<MyTraits>;

        std::unique_ptr<SimulationController> sim;

        static constexpr auto argv = std::array{"hemelb", "-in", "large_cylinder_rbc.xml"};
        configuration::CommandLine options = {argv.size(), argv.data()};
        redblood::VTKMeshIO vtk_io = {};

    public:
        MulticellBenchmarkTests() : FolderTestFixture() {
            CopyResourceToTempdir("large_cylinder_rbc.xml");
            CopyResourceToTempdir("large_cylinder.gmy");
            CopyResourceToTempdir("rbc_ico_720.msh");

            ModifyXMLInput("large_cylinder_rbc.xml", { "simulation", "steps", "value" }, 64000); // Enough time for 6 cells to coexist in the domain
            ModifyXMLInput("large_cylinder_rbc.xml",
                           { "redbloodcells", "controller", "stencil" },
                           "two");
            ModifyXMLInput("large_cylinder_rbc.xml", { "redbloodcells",
                                                       "templates",
                                                       "cell",
                                                       "shape",
                                                       "mesh_path" },
                           "rbc_ico_720.msh");
            ModifyXMLInput("large_cylinder_rbc.xml", { "inlets",
                                                       "inlet",
                                                       "insertcell",
                                                       "every",
                                                       "value" },
                           2.5e-3); // 2.5e-3 is 12500 timesteps
            ModifyXMLInput("large_cylinder_rbc.xml", { "inlets",
                                                       "inlet",
                                                       "insertcell",
                                                       "delta_x",
                                                       "value" },
                           3e-6);
            ModifyXMLInput("large_cylinder_rbc.xml", { "inlets",
                                                       "inlet",
                                                       "insertcell",
                                                       "delta_y",
                                                       "value" },
                           3e-6);

            sim = configuration::SimBuilder::CreateSim<MyTraits>(options, Comms());
        }

        void testBenchmark()
        {
            unsigned timestep = 0;
            auto output_callback =
                    [this, &timestep](const redblood::CellContainer & cells) {
                        if ((timestep % 100) == 0) {
                            for (auto& cell: cells) {
                                std::stringstream filename;
                                filename << cell->GetTag() << "_t_" << timestep << ".vtp";
                                vtk_io.writeFile(filename.str(), *cell, &this->sim->GetUnitConverter());
                            }
                        }
                        timestep++;
                    };

#ifdef HEMELB_CALLGRIND
            bool instrumented = false;
            unsigned int start_timestep;
            auto callgrind_callback =
                    [&timestep, &instrumented, &start_timestep](const hemelb::redblood::CellContainer & cells) {
                        if (cells.size() >= 5 && !instrumented) {
                            instrumented = true;
                            start_timestep = timestep;
                            std::cerr << "Starting callgrind instrumentation at timestep " << timestep << std::endl;
                            CALLGRIND_ZERO_STATS;
                            CALLGRIND_START_INSTRUMENTATION;
                        }
                        if (instrumented && timestep - start_timestep == 1000) {
                            CALLGRIND_STOP_INSTRUMENTATION;
                            std::cerr << "Stopped callgrind instrumentation at timestep " << timestep << std::endl;
                            CALLGRIND_DUMP_STATS;
                        }
                    };
#endif

            REQUIRE(sim);
            auto controller = std::static_pointer_cast<CellControl>(sim->GetCellController());
            REQUIRE(controller);
            controller->AddCellChangeListener(output_callback);
#ifdef HEMELB_CALLGRIND
            controller->AddCellChangeListener(callgrind_callback);
#endif

            // run the simulation
            sim->RunSimulation();
            sim->Finalise();

            AssertPresent("results/report.txt");
            AssertPresent("results/report.xml");
        }

    };

#ifdef HEMELB_CALLGRIND
    METHOD_AS_TEST_CASE(MulticellBenchmarkTests::testBenchmark,
			"Run a simulation long enought to have 6 cells",
			"[redblood][.long]");
#endif
}
