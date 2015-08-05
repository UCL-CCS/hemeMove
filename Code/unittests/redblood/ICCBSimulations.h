#ifndef HEMELB_UNITTESTS_REDBLOOD_ICCBSIMULATIONS_H
#define HEMELB_UNITTESTS_REDBLOOD_ICCBSIMULATIONS_H

#include <cppunit/extensions/HelperMacros.h>
#include <memory>

#include "lb/BuildSystemInterface.h"
#include "Traits.h"
#include "redblood/Mesh.h"
#include "redblood/Cell.h"
#include "redblood/CellController.h"
#include "unittests/helpers/FolderTestFixture.h"
#include <boost/uuid/uuid_io.hpp>

namespace hemelb
{
  namespace unittests
  {
    namespace redblood
    {
      class ICCBSimulations : public hemelb::unittests::helpers::FolderTestFixture
      {
          CPPUNIT_TEST_SUITE (ICCBSimulations);
          CPPUNIT_TEST (testSimpleTube);
          CPPUNIT_TEST_SUITE_END();

          typedef Traits<>::ChangeKernel<lb::GuoForcingLBGK>::Type Traits;
          typedef hemelb::redblood::CellController<Traits::Kernel> CellControl;
          typedef SimulationMaster<Traits> MasterSim;

        public:
          void setUp()
          {
            FolderTestFixture::setUp();
            CopyResourceToTempdir("large_cylinder_rbc.xml");
            CopyResourceToTempdir("large_cylinder.gmy");
            CopyResourceToTempdir("red_blood_cell.txt");

            std::vector<std::string> intel;
            intel.push_back("simulation");
            intel.push_back("steps");
            intel.push_back("value");
            ModifyXMLInput("large_cylinder_rbc.xml", std::move(intel), 300000);

            argv[0] = "hemelb";
            argv[1] = "-in";
            argv[2] = "large_cylinder_rbc.xml";
            argv[3] = "-i";
            argv[4] = "1";
            argv[5] = "-ss";
            argv[6] = "1111";
            options = std::make_shared<hemelb::configuration::CommandLine>(argc, argv);

            master = std::make_shared<MasterSim>(*options, Comms());
          }

          void tearDown()
          {
            master->Finalise();
            master.reset();
          }

          void testSimpleTube()
          {
            unsigned timestep = 0;
            auto output_callback = [this, &timestep](const hemelb::redblood::CellContainer & cells)
            {
              for (auto cell: cells)
              {
                std::stringstream filename;
                filename << cell->GetTag() << "_t_" << timestep++ << ".vtp";
                hemelb::redblood::writeVTKMesh(filename.str(), cell, this->master->GetUnitConverter());
              }

            };
            CPPUNIT_ASSERT(master);
            auto controller = std::static_pointer_cast<CellControl>(master->GetCellController());
            CPPUNIT_ASSERT(controller);
            controller->AddCellChangeListener(output_callback);

            // run the simulation
            master->RunSimulation();

            AssertPresent("results/report.txt");
            AssertPresent("results/report.xml");
          }

        private:
          std::shared_ptr<MasterSim> master;
          std::shared_ptr<hemelb::configuration::CommandLine> options;
          int const argc = 7;
          char const * argv[7];

      };
      // Extra lines needed by CPPUNIT_TEST_SUITE_REGISTRATION
      // Extra lines needed by CPPUNIT_TEST_SUITE_REGISTRATION
      // Don't register the unit test so it's not run by default as part of CI.
      // Uncomment the line below in order to run the test with:
      //   ./unittests_hemelb hemelb::unittests::redblood::ICCBSimulations::testSimpleTube
      //CPPUNIT_TEST_SUITE_REGISTRATION (ICCBSimulations);

    } // namespace redblood
  } // namespace unittests
} // namespace hemelb

#endif
