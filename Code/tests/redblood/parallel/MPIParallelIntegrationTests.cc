// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <random>
#include <memory>

#include <catch2/catch.hpp>

#include "redblood/CellController.h"
#include "redblood/parallel/IntegrateVelocities.h"
#include "redblood/parallel/CellParallelization.h"
#include "configuration/CommandLine.h"
#include "tests/redblood/Fixtures.h"
#include "tests/helpers/LatticeDataAccess.h"
#include "tests/helpers/FolderTestFixture.h"
#include "tests/helpers/ApproxVector.h"
#include "tests/redblood/parallel/ParallelFixtures.h"

namespace hemelb::tests
{
    using namespace redblood;

    class MPIParallelIntegrationTests : public OpenSimFixture
    {
    public:
        MPIParallelIntegrationTests();

        template<class STENCIL> void testIntegration()
        {
            Check<STENCIL>();
        }

    protected:
        template<class STENCIL> void Check();
    };

    MPIParallelIntegrationTests::MPIParallelIntegrationTests() : OpenSimFixture()
    {
      // Have everything ready to creates simulations
      if (net::MpiCommunicator::World().Rank() == 0) {
	CopyResourceToTempdir("large_cylinder_rbc.xml");
	CopyResourceToTempdir("large_cylinder.gmy");
	CopyResourceToTempdir("red_blood_cell.txt");

	// Get the RBC drop off point closer to a subdomain border so
	// that we can test for cell communication related issues with
	// fewer timesteps.
	ModifyXMLInput("large_cylinder_rbc.xml",
		       { "simulation", "steps", "value" },
		       5000);
	ModifyXMLInput("large_cylinder_rbc.xml",
		       { "inlets", "inlet", "flowextension", "length", "value" },
		       20e-6);

      }
      net::MpiCommunicator::World().Barrier();

      auto const result = Comms().Rank() == 0 ?
	"root_result" :
	"others";
      options = std::make_shared<configuration::CommandLine>(configuration::CommandLine { "hemelb",
	    "-in",
	    "large_cylinder_rbc.xml",
	    "-out",
	    result });
    }


    template<class STENCIL>
    void MPIParallelIntegrationTests::Check() {
      using hemelb::redblood::CellContainer;
      using hemelb::redblood::TemplateCellContainer;
      using T = MyTraits<STENCIL>;

      auto const world = net::MpiCommunicator::World();
      auto const color = world.Rank() == 0;
      auto const split = net::IOCommunicator(world.Split(color));

      auto sim = CreateSim<STENCIL>(split);
      REQUIRE(sim);

      unsigned num_cells;
      auto checkNumCells = [&num_cells]( const CellContainer & cells) {
	num_cells = cells.size();
      };

      auto controller = std::static_pointer_cast<CellController<T>>(sim->GetCellController());
      REQUIRE(controller);
      controller->AddCellChangeListener(checkNumCells);

      // run the simulation
      sim->RunSimulation();
      sim->Finalise();

      // check that both simulations have one cell in the domain
      unsigned num_cells_sequential;
      unsigned num_cells_parallel;
      if (color) {
	num_cells_sequential = num_cells;
	num_cells_parallel = 0;
      } else {
	num_cells_parallel = num_cells;
      }
      world.Broadcast(num_cells_sequential, 0);
      num_cells_parallel = world.AllReduce(num_cells_parallel, MPI_SUM);
      REQUIRE(num_cells_sequential == 1u);
      REQUIRE(num_cells_parallel == 1u);

      // Find out number of vertices in the RBC mesh
      unsigned num_vertices = color ? (*controller->GetCells().begin())->GetVertices().size() : 0;
      world.Broadcast(num_vertices, 0);

      // Find out the location of the mesh vertices in the sequential run
      hemelb::redblood::MeshData::Vertices sequential_vertices(num_vertices);
      if (color) {
	sequential_vertices = (*controller->GetCells().begin())->GetVertices();
      }
      HEMELB_MPI_CALL(MPI_Bcast, (&sequential_vertices[0], 3*sequential_vertices.size(), net::MpiDataType<LatticeDistance>(), 0, net::MpiCommunicator::World()));

      // Find out the location of the mesh vertices in the parallel run
      hemelb::redblood::MeshData::Vertices vertices_location_reduction(num_vertices, {0,0,0});
      if (!color && (controller->GetCells().size()==1)) {
	// If I'm the process in the parallel run than owns the only existing cell
	vertices_location_reduction = (*controller->GetCells().begin())->GetVertices();
      }
      hemelb::redblood::MeshData::Vertices parallel_vertices(num_vertices);
      HEMELB_MPI_CALL(MPI_Allreduce,
		      (vertices_location_reduction.data(), parallel_vertices.data(), 3*vertices_location_reduction.size(), net::MpiDataType<LatticeDistance>(), MPI_SUM, net::MpiCommunicator::World()));

      // Compare locations in sequential and parallel runs
      auto approx = ApproxVector<LatticePosition>(0.0).Margin(1e-12);
      for (auto const item : util::zip(sequential_vertices, parallel_vertices)) {
	REQUIRE(approx(std::get<0>(item)) == std::get<1>(item));
	// auto delta = 1e-12;
	// REQUIRE_DOUBLES_EQUAL(std::get<0>(item)[0], std::get<1>(item)[0], delta);
	// REQUIRE_DOUBLES_EQUAL(std::get<0>(item)[1], std::get<1>(item)[1], delta);
	// REQUIRE_DOUBLES_EQUAL(std::get<0>(item)[2], std::get<1>(item)[2], delta);
      }

    }
    METHOD_AS_TEST_CASE(MPIParallelIntegrationTests::testIntegration<hemelb::redblood::stencil::FourPoint>,
			"MPIParallelIntegrationTests with FourPoint stencil",
			"[redblood][.long]");
  }
