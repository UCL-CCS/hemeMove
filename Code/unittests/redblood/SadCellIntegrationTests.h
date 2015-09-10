#ifndef HEMELB_UNITTESTS_REDBLOOD_SADCELL_INTEGRATION_H
#define HEMELB_UNITTESTS_REDBLOOD_SADCELL_INTEGRATION_H

#include <cppunit/extensions/HelperMacros.h>
#include <boost/uuid/uuid_io.hpp>
#include <memory>

#include "SimulationMaster.h"
#include "lb/BuildSystemInterface.h"
#include "lb/lattices/D3Q19.h"
#include "Traits.h"
#include "redblood/Mesh.h"
#include "redblood/Cell.h"
#include "redblood/CellController.h"
#include "unittests/helpers/FolderTestFixture.h"

namespace hemelb
{
  namespace unittests
  {
    namespace redblood
    {
      class SadCellIntegrationTests : public hemelb::unittests::helpers::FolderTestFixture
      {
          CPPUNIT_TEST_SUITE (SadCellIntegrationTests);
          CPPUNIT_TEST (testIntegration);CPPUNIT_TEST_SUITE_END()
          ;

          typedef Traits<>::Reinstantiate<lb::lattices::D3Q19, lb::GuoForcingLBGK>::Type Traits;
          typedef hemelb::redblood::CellController<Traits::Kernel> CellControl;
          typedef SimulationMaster<Traits> MasterSim;

        public:
          void setUp()
          {
            FolderTestFixture::setUp();
            CopyResourceToTempdir("large_cylinder_rbc.xml");
            CopyResourceToTempdir("large_cylinder.gmy");
            CopyResourceToTempdir("rbc_ico_1280.msh");

            ModifyXMLInput(
                "large_cylinder_rbc.xml", {"simulation", "steps", "value"}, 50000);
            ModifyXMLInput(
                "large_cylinder_rbc.xml", {"redbloodcells", "controller", "stencil"}, "two");
            ModifyXMLInput(
                "large_cylinder_rbc.xml", {"redbloodcells", "cells", "cell", "shape", "mesh_path"},
                "rbc_ico_1280.msh");
            ModifyXMLInput(
                "large_cylinder_rbc.xml", {"inlet", "condition", "mean", "value"}, 0);
            DeleteXMLInput("large_cylinder_rbc.xml", {"redbloodcells", "cells"});
            DeleteXMLInput("large_cylinder_rbc.xml", {"inlets", "inlet", "flowextension"});
            DeleteXMLInput("large_cylinder_rbc.xml", {"outlets", "outlet", "flowextension"});

            argv[0] = "hemelb";
            argv[1] = "-in";
            argv[2] = "large_cylinder_rbc.xml";
            argv[3] = "-i";
            argv[4] = "0";
            argv[5] = "-ss";
            argv[6] = "1111";
            options = std::make_shared<hemelb::configuration::CommandLine>(argc, argv);

            master = std::make_shared<MasterSim>(*options, Comms());
          }

          void testIntegration()
          {
            auto const normal = readMesh(resources::Resource("rbc_ico_1280.msh").Path().c_str());
            auto const cell = std::make_shared<redblood::Cell>(normal->vertices, normal);
            auto const & converter = master->GetUnitConverter();
            auto const deformed = readMesh(resources::Resource("sad.msh").Path().c_str());
            auto const sadcell = std::make_shared<redblood::Cell>(deformed->vertices, normal);
            auto const scale = converter.ConvertToLatticeUnits("m", 4e-6);
            // scale and positions cell in the middle somewhere
            cell->SetScale(scale);
            *cell *= scale;
            sadcell->SetScale(scale);
            *sadcell *= 1e0/converter.GetVoxelSize();
            std::cout << sadcell->GetVertices()[0] << "\n"
                      << sadcell->GetVertices()[1] << "\n"
                      << sadcell->GetVertices()[2] << "\n";
            *sadcell +=
              converter.ConvertPositionToLatticeUnits(PhysicalPosition(0, 0, 0))
              - sadcell->GetBarycenter();
            *cell +=
              converter.ConvertPositionToLatticeUnits(PhysicalPosition(0, 0, 0))
              - cell->GetBarycenter();
            writeVTKMesh("/tmp/ideal.vtp", cell, converter);
            writeVTKMesh("/tmp/deformed.vtp", sadcell, converter);
            sadcell->moduli.bending = converter.ConvertToLatticeUnits("Nm", 2e-19);
            sadcell->moduli.surface = 1e0;
            sadcell->moduli.volume = 1e0;
            sadcell->moduli.dilation = 0.75;
            sadcell->moduli.strain = converter.ConvertToLatticeUnits("N/m", 5e-6);

            std::cout << "VOLUMES: " << sadcell->GetVolume() << "  " << cell->GetVolume() << "\n";
            auto controller = std::static_pointer_cast<CellControl>(master->GetCellController());
            controller->AddCell(sadcell);
            std::vector<PhysicalEnergy> energies;
            energies.push_back((*sadcell)());
            controller->AddCellChangeListener(
                [&energies](CellContainer const &cells)
                {
                   energies.push_back((**cells.begin())());
                }
            );
            // run the simulation
            master->RunSimulation();
            for(auto const energy: energies)
            {
              std::cout << std::scientific << " energy: " << energy << std::endl;
            }
            writeVTKMesh("/tmp/reformed.vtp", sadcell, converter);

            *sadcell +=
              converter.ConvertPositionToLatticeUnits(PhysicalPosition(0, 0, 0))
              - sadcell->GetBarycenter();
            writeVTKMesh("/tmp/reformed_centered.vtp", sadcell, converter);

            AssertPresent("results/report.txt");
            AssertPresent("results/report.xml");
          }

        private:
          std::shared_ptr<MasterSim> master;
          std::shared_ptr<hemelb::configuration::CommandLine> options;
          int const argc = 7;
          char const * argv[7];

      };






      CPPUNIT_TEST_SUITE_REGISTRATION (SadCellIntegrationTests);
    } // namespace redblood
  } // namespace unittests
} // namespace hemelb

#endif
