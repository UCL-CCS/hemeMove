// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <catch2/catch.hpp>
#include <tinyxml2.h>

#include <algorithm>
#include <random>

#include "redblood/parallel/SpreadForces.h"
#include "configuration/CommandLine.h"
#include "SimulationController.h"
#include "util/span.h"
#include "tests/redblood/Fixtures.h"
#include "tests/helpers/LatticeDataAccess.h"
#include "tests/helpers/FolderTestFixture.h"
#include "tests/redblood/parallel/ParallelFixtures.h"
#include "tests/helpers/ApproxVector.h"

namespace hemelb::tests
{
    class MPISpreadForcesTests : public OpenSimFixture
    {
    public:
        MPISpreadForcesTests();

        template<class STENCIL> void testMidRegion()
        {
            Check<STENCIL>(2, 0, 1);
        }
        template<class STENCIL> void testEdgeRegion()
        {
            Check<STENCIL>(0, 2, 1);
        }
        template<class STENCIL> void testAll()
        {
            Check<STENCIL>(10, 10, 2);
        }

    protected:
        template<class STENCIL> void Check(size_t mid, size_t edges, size_t nCells);
    };

    MPISpreadForcesTests::MPISpreadForcesTests() : OpenSimFixture()
    {
      using hemelb::configuration::CommandLine;

      // Have everything ready to creates simulations
      if (net::MpiCommunicator::World().Rank() == 0) {
	CopyResourceToTempdir("red_blood_cell.txt");
	tinyxml2::XMLDocument doc;
    doc.LoadFile(resources::Resource("large_cylinder.xml").Path().c_str());
	CopyResourceToTempdir("large_cylinder.xml");
	ModifyXMLInput("large_cylinder.xml", { "simulation", "steps", "value" }, 2);
	CopyResourceToTempdir("large_cylinder.gmy");
      }
      net::MpiCommunicator::World().Barrier();
      auto const result = Comms().Rank() == 0 ?
	"root_result" :
	"others";
      options = std::make_shared<CommandLine>(CommandLine { "hemelb",
	    "-in",
	    "large_cylinder.xml",
	    "-out",
	    result });
    }

    template<class STENCIL>
    void MPISpreadForcesTests::Check(size_t mid, size_t edges, size_t nCells)
    {
      using hemelb::redblood::CellContainer;
      using hemelb::redblood::parallel::nodeDistributions;
      auto const world = net::MpiCommunicator::World();
      if (world.Size() == 1)
        {
          return;
        }
      auto const color = world.Rank() == 0;
      auto const split = net::IOCommunicator(world.Split(color));
      auto sim = CreateSim<STENCIL>(split);
      auto& fieldData = sim->GetFieldData();
      auto& dom = fieldData.GetDomain();
      helpers::ZeroOutForces(fieldData);

      // Figure out positions to use for cell nodes
      auto const cells = CreateCellsFromSpecialPositions(dom, mid, edges, world, nCells);
      auto const owned =
	split.Size() != 1 ?
	CellContainer { cells.begin() + split.Rank() * nCells, cells.begin()
			+ (1 + split.Rank()) * nCells } :
      CellContainer { cells.begin(), cells.end() };
      auto const& graphComm = CreateDumbGraphComm(split);
      auto const distributions = nodeDistributions(hemelb::redblood::parallel::ComputeGlobalCoordsToProcMap(graphComm, dom), owned);

      hemelb::redblood::parallel::SpreadForces mpi_spreader(graphComm);
      mpi_spreader.PostMessageLength(distributions, owned);
      mpi_spreader.ComputeForces(owned);
      mpi_spreader.PostForcesAndNodes(distributions, owned);
      mpi_spreader.SpreadLocalForces<MyTraits<STENCIL>>(fieldData, owned);
      mpi_spreader.SpreadNonLocalForces<MyTraits<STENCIL>>(fieldData);

      std::vector<LatticeVector> indices;
      std::vector<LatticeForceVector> forces;
      if (color)
        {
          for (site_t i = 0; i < dom.GetLocalFluidSiteCount(); ++i)
	    {
	      auto const site = fieldData.GetSite(i);
	      if (site.GetForce().GetMagnitudeSquared() > 1e-8)
		{
		  indices.push_back(site.GetGlobalSiteCoords());
		  forces.push_back(site.GetForce());
		}
	    }
        }
      int nIndices = indices.size();
      world.Broadcast(nIndices, 0);
      REQUIRE(nIndices > 0);
      indices.resize(nIndices);
      world.Broadcast(to_span(indices), 0);

      if (not color)
        {
          for (auto coords : indices)
	    {
	      auto const id = dom.GetProcIdFromGlobalCoords(coords);
	      if (id == split.Rank())
		{
		  auto const site = fieldData.GetSite(coords);
		  forces.push_back(-site.GetForce());
		}
	      else
		{
		  forces.push_back(LatticeForceVector::Zero());
		}
	    }
        }
      REQUIRE(size_t(nIndices) == size_t(forces.size()));
      std::vector<LatticeForceVector> summed(indices.size());
      HEMELB_MPI_CALL(MPI_Allreduce,
		      (forces.data(), summed.data(), forces.size() *
        3, net::MpiDataType<double>(), MPI_SUM, world));
      auto zero = ApproxVector<LatticeForceVector>(0.0).Margin(1e-8);
      for (auto const& force : summed) {
	REQUIRE(force == zero);
      }
    }

    METHOD_AS_TEST_CASE(MPISpreadForcesTests::testMidRegion<hemelb::redblood::stencil::FourPoint>,
			"Test mid region with 4 point stencil",
			"[redblood]");
    METHOD_AS_TEST_CASE(MPISpreadForcesTests::testEdgeRegion<hemelb::redblood::stencil::FourPoint>,
			"Test edge region with 4 point stencil",
			"[redblood]");
    METHOD_AS_TEST_CASE(MPISpreadForcesTests::testMidRegion<hemelb::redblood::stencil::FourPoint>,
			"Test all with 4 point stencil",
			"[redblood]");

    METHOD_AS_TEST_CASE(MPISpreadForcesTests::testMidRegion<hemelb::redblood::stencil::ThreePoint>,
			"Test mid region with 3 point stencil",
			"[redblood]");
    METHOD_AS_TEST_CASE(MPISpreadForcesTests::testEdgeRegion<hemelb::redblood::stencil::ThreePoint>,
			"Test edge region with 3 point stencil",
			"[redblood]");
    METHOD_AS_TEST_CASE(MPISpreadForcesTests::testMidRegion<hemelb::redblood::stencil::ThreePoint>,
			"Test all with 3 point stencil",
			"[redblood]");

    METHOD_AS_TEST_CASE(MPISpreadForcesTests::testMidRegion<hemelb::redblood::stencil::CosineApprox>,
			"Test mid region with cosine stencil",
			"[redblood]");
    METHOD_AS_TEST_CASE(MPISpreadForcesTests::testEdgeRegion<hemelb::redblood::stencil::CosineApprox>,
			"Test edge region with cosine stencil",
			"[redblood]");
    METHOD_AS_TEST_CASE(MPISpreadForcesTests::testMidRegion<hemelb::redblood::stencil::CosineApprox>,
			"Test all with cosine stencil",
			"[redblood]");

    METHOD_AS_TEST_CASE(MPISpreadForcesTests::testMidRegion<hemelb::redblood::stencil::TwoPoint>,
			"Test mid region with 2 point stencil",
			"[redblood]");
    METHOD_AS_TEST_CASE(MPISpreadForcesTests::testEdgeRegion<hemelb::redblood::stencil::TwoPoint>,
			"Test edge region with 2 point stencil",
			"[redblood]");
    METHOD_AS_TEST_CASE(MPISpreadForcesTests::testMidRegion<hemelb::redblood::stencil::TwoPoint>,
			"Test all with 2 point stencil",
			"[redblood]");

  }

