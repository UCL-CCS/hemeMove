// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <random>
#include <memory>

#include <catch2/catch.hpp>

#include "util/span.h"
#include "redblood/parallel/NodeCharacterizer.h"
#include "configuration/CommandLine.h"
#include "configuration/SimBuilder.h"
#include "tests/redblood/Fixtures.h"
#include "tests/helpers/LatticeDataAccess.h"
#include "tests/helpers/FolderTestFixture.h"
#include "tests/redblood/parallel/ParallelFixtures.h"

namespace hemelb::tests
{
    using namespace redblood;
    //! Parallel and Sequential move in lock step
    class ParallelFixtureTests : public OpenSimFixture
    {
    public:
        ParallelFixtureTests();

        //! if owner procs thinks position affects proc i, then proc i knows it as well
        void testTransitiveOwnership();
    };


    ParallelFixtureTests::ParallelFixtureTests() : OpenSimFixture()
    {
      // Have everything ready to create simulations
      if (Comms().Rank() == 0)
	{
	  CopyResourceToTempdir("large_cylinder_rbc.xml");
	  CopyResourceToTempdir("large_cylinder.gmy");
	  CopyResourceToTempdir("red_blood_cell.txt");

	  ModifyXMLInput("large_cylinder_rbc.xml", { "simulation", "steps", "value" }, 20000);
	}
      Comms().Barrier();

      options = std::make_shared<configuration::CommandLine>(configuration::CommandLine { "hemelb",
	    "-in",
	    "large_cylinder_rbc.xml",
	    "-out",
	    "" });
    }

    void ParallelFixtureTests::testTransitiveOwnership()
    {
      using Stencil = stencil::FourPoint;

      auto const world = Comms();
      if(world.Size() < 2)
        {
          return;
        }

      auto sim = CreateSim<Stencil>(world);
      REQUIRE(sim);

      auto& fd = sim->GetFieldData();
      auto& dom = fd.GetDomain();
      helpers::ZeroOutForces(fd);
      auto const nmid = 20;
      auto const nedges = 20;
      auto const positions = GatherSpecialPositions(dom, nmid, nedges, world);

      auto graphComm =
	world.DistGraphAdjacent(parallel::ComputeProcessorNeighbourhood(world,
                                                                          dom,
									      2e-6 / sim->GetSimConfig().GetVoxelSize()));

      auto const& globalCoordsToProcMap = parallel::ComputeGlobalCoordsToProcMap(graphComm, dom);

      for(std::size_t i(0); i < positions.size(); ++i)
        {
          auto const procs = parallel::details::positionAffectsProcs<Stencil>(
												globalCoordsToProcMap, positions[i].as<LatticeDistance>());

          // Send set of affected procs as known by owner proc
          decltype(world.Rank()) positions_are_from_this_proc = i / (nmid + nedges) + 1;
          int N(procs.size());
          world.Broadcast(N, positions_are_from_this_proc);
          std::vector<proc_t> expected(N);
          if(world.Rank() == positions_are_from_this_proc)
	    {
	      std::copy(procs.begin(), procs.end(), expected.begin());
	    }
          world.Broadcast(to_span(expected), positions_are_from_this_proc);

          std::set<proc_t> expected_set(expected.begin(), expected.end());

          // Owner knows thinks current position affects this proc
          if(expected_set.count(world.Rank()))
	    {
	      // so this proc must know that current positions affects it
	      REQUIRE(procs.count(world.Rank()));
	      // and affects owner proc
	      REQUIRE(procs.count(positions_are_from_this_proc));
	    }
          else
	    {
	      // otherwise, this proc should not think it is affected
	      REQUIRE(not procs.count(world.Rank()));
	    }

          // similarly, if this proc thinks current position affects owner proc
          if(not procs.count(positions_are_from_this_proc))
	    {
	      // then owner proc must think current position affects this proc
	      REQUIRE(not expected_set.count(world.Rank()));
	    }

        }
    }

    METHOD_AS_TEST_CASE(ParallelFixtureTests::testTransitiveOwnership,
			"Test transitive ownership",
			"[redblood]");
}

