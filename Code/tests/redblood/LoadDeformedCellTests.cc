// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <memory>

#include <catch2/catch.hpp>
#include <boost/uuid/uuid_io.hpp>

#include "configuration/SimBuilder.h"
#include "lb/lattices/D3Q19.h"
#include "Traits.h"
#include "redblood/Mesh.h"
#include "redblood/Cell.h"
#include "redblood/CellController.h"
#include "redblood/MeshIO.h"

#include "tests/helpers/FolderTestFixture.h"

namespace hemelb::tests
  {
    using namespace redblood;
    TEST_CASE_METHOD(helpers::FolderTestFixture, "LoadDeformedCellTests") {

      using Traits = Traits<lb::D3Q19, lb::GuoForcingLBGK>;
      using CellControl = hemelb::redblood::CellController<Traits>;

      redblood::VTKMeshIO io = {};

      CopyResourceToTempdir("large_cylinder_rbc.xml");
      CopyResourceToTempdir("large_cylinder.gmy");
      CopyResourceToTempdir("rbc_ico_720.vtp");
      CopyResourceToTempdir("992Particles_rank3_26_t992.vtp");

      // Run simulation for longer with no flow across the cylinder
      ModifyXMLInput("large_cylinder_rbc.xml", { "simulation", "steps", "value" }, 1000);
      ModifyXMLInput("large_cylinder_rbc.xml", { "inlets",
            "inlet",
            "condition",
            "mean",
            "value" },
        0);

      // Remove cell insertion/removal and configuration
      DeleteXMLInput("large_cylinder_rbc.xml", { "inlets", "inlet", "insertcell" });
      DeleteXMLInput("large_cylinder_rbc.xml", { "inlets", "inlet", "flowextension" });
      DeleteXMLInput("large_cylinder_rbc.xml", { "outlets", "outlet", "flowextension" });
      DeleteXMLInput("large_cylinder_rbc.xml", { "redbloodcells", "templates" });

      int constexpr argc = 3;
      char const * argv[argc] = {
        "hemelb", "-in", "large_cylinder_rbc.xml",
      };
      auto options = std::make_shared<configuration::CommandLine>(argc, argv);
      auto sim = configuration::SimBuilder::CreateSim<Traits>(*options, Comms());

      SECTION("testIntegration") {
        // Read meshes from disc
        auto const normal = io.readFile("rbc_ico_720.vtp", true);
        auto const deformed = io.readFile("992Particles_rank3_26_t992.vtp", true);

        // Check they are compatible
        REQUIRE(normal->facets == deformed->facets);
        REQUIRE(volume(*normal) > 0.0);
        REQUIRE(volume(deformed->vertices, normal->facets) > 0.0);

        auto const& converter = sim->GetUnitConverter();
        auto const scale = converter.ConvertToLatticeUnits("m", 4e-6);
        auto const cell = std::make_shared<redblood::Cell>(normal->vertices,
                                                           Mesh(normal),
                                                           scale);
        auto const sadcell = std::make_shared<redblood::Cell>(deformed->vertices,
                                                              Mesh(normal),
                                                              scale);

        // Scale and move to the centre of the domain
        *cell *= scale;
        *sadcell *= scale;
        *cell += converter.ConvertPositionToLatticeUnits(PhysicalPosition(0, 0, 0))
          - cell->GetBarycentre();
        *sadcell += converter.ConvertPositionToLatticeUnits(PhysicalPosition(0, 0, 0))
          - sadcell->GetBarycentre();

        io.writeFile("ideal.vtp", *cell, &converter);
        io.writeFile("deformed.vtp", *sadcell, &converter);

        sadcell->moduli.bending = 0.1;
        sadcell->moduli.strain = 0.1;
        sadcell->moduli.surface = 1e0;
        sadcell->moduli.volume = 1e0;
        sadcell->moduli.dilation = 0.75;

        auto controller = std::static_pointer_cast<CellControl>(sim->GetCellController());
        controller->AddCell(sadcell);

        controller->AddCellChangeListener(
          [&converter,io](const hemelb::redblood::CellContainer &cells)
          {
            static int iter = 0;
            for (auto& cell : cells) {
              if(iter % 10 == 0) {
                std::stringstream filename;
                filename << cell->GetTag() << "_t_" << iter << ".vtp";
                io.writeFile(filename.str(), *cell, &converter);
              }
            }
            ++iter;
          });

        // run the simulation
        sim->RunSimulation();

        // Check simulation ran until the end
        AssertPresent("results/report.txt");
        AssertPresent("results/report.xml");

        // Recentre simulated cell
        io.writeFile("reformed.vtp", *sadcell, &converter);
        *sadcell += converter.ConvertPositionToLatticeUnits(PhysicalPosition(0, 0, 0))
          - sadcell->GetBarycentre();
        io.writeFile("reformed_centered.vtp", *sadcell, &converter);

//            // TODO: Align both cells for comparison
//            auto cell01 = cell->GetVertices()[350] - cell->GetVertices()[154];
//            auto sadcell01 = sadcell->GetVertices()[350] - sadcell->GetVertices()[154];
//            auto rotation = rotationMatrix(sadcell01, cell01);
//            *sadcell *= rotation;
//            writeVTKMesh("/tmp/reformed_centered_reorientied.vtp", sadcell, converter);
//            writeVTKMesh("/tmp/cell_end.vtp", cell, converter);
//
//            auto compareVertices = [] (util::Vector3D<double> const &left, util::Vector3D<double> const &right){
//              std::cout << left << std::endl;
//              std::cout << right << std::endl << std::endl;
//              return (left-right).GetMagnitudeSquared() < 1e-4;};
//            CPPUNIT_ASSERT(std::equal(cell->GetVertices().begin(), cell->GetVertices().end(), sadcell->GetVertices().begin(), compareVertices));
      }
    }

  } // namespace hemelb
