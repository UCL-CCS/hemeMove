// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include<algorithm>
#include <random>

#include <catch2/catch.hpp>
#include <tinyxml2.h>

#include "redblood/parallel/IntegrateVelocities.h"
#include "redblood/parallel/CellParallelization.h"
#include "configuration/CommandLine.h"
#include "SimulationController.h"
#include "util/span.h"
#include "tests/redblood/Fixtures.h"
#include "tests/helpers/LatticeDataAccess.h"
#include "tests/helpers/FolderTestFixture.h"
#include "tests/helpers/ApproxVector.h"
#include "tests/redblood/parallel/ParallelFixtures.h"

namespace hemelb::tests
{
    using namespace redblood;

    class MPIIntegrateVelocitiesTests : public OpenSimFixture
    {
    public:
        MPIIntegrateVelocitiesTests();

        template<class STENCIL>
        void testMidRegion()
        {
            Check<STENCIL>(2, 0, 1);
        }
        template<class STENCIL>
        void testEdgeRegion()
        {
            Check<STENCIL>(0, 2, 1);
        }
        template<class STENCIL>
        void testAll()
        {
            Check<STENCIL>(3, 3, 2);
        }

    protected:
        template<class STENCIL>
        void Check(size_t mid, size_t edges, size_t nCells);
    };

    MPIIntegrateVelocitiesTests::MPIIntegrateVelocitiesTests() : OpenSimFixture()
    {
      using hemelb::configuration::CommandLine;

      // Have everything ready to creates simulations
      if (net::MpiCommunicator::World().Rank() == 0)
        {
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
    void MPIIntegrateVelocitiesTests::Check(size_t mid, size_t edges, size_t nCells)
    {
      using hemelb::redblood::CellContainer;
      using hemelb::redblood::TemplateCellContainer;
      using Traits = MyTraits<STENCIL>;

      auto const world = net::MpiCommunicator::World();
      auto const color = world.Rank() == 0;
      auto const split = net::IOCommunicator(world.Split(color));
      auto sim = CreateSim<STENCIL>(split);
      auto& fieldData = sim->GetFieldData();
      auto& dom = fieldData.GetDomain();
      helpers::ZeroOutForces(fieldData);

      // Setup same distribution for all instances
      for (Direction i = 0; i < Traits::Lattice::NUMVECTORS; ++i) {
	auto func = [i](LatticePosition const &x) {
	  double const fac = (1e0 + static_cast<double>(i) * 0.1);
	  return x.GetMagnitudeSquared() * fac - Dot(x, LatticePosition{1, 1, 1});
	};
	helpers::setUpDistribution<typename Traits::Lattice>(&fieldData, i, func);
      }

      // Figure out positions to use for cell nodes
      auto const cells = CreateCellsFromSpecialPositions(dom, mid, edges, world, nCells);
      auto const checkMoves = cells[0]->GetVertices()[0];
      auto owned =
	split.Size() != 1 ?
	CellContainer { cells.begin() + split.Rank() * nCells, cells.begin()
			+ (1 + split.Rank()) * nCells } :
      CellContainer { cells.begin(), cells.end() };
      auto const& graphComm = CreateDumbGraphComm(split);
      auto const distributions = hemelb::redblood::parallel::nodeDistributions(hemelb::redblood::parallel::ComputeGlobalCoordsToProcMap(graphComm, dom), owned);

      // Goes through "ExchangeCells" to figure out who owns/lends what.
      // Ownership is pre-determined here: first nCells got to 0, second nCells to 2, etc...
      TemplateCellContainer const templates =
	{ { cells[0]->GetTemplateName(), cells[0]->clone() } };
      auto ownership = [&cells, nCells](CellContainer::const_reference cell) {
	size_t i(0);
	for(auto const& c: cells) {
	  if(cell->GetTag() == c->GetTag()) {
	    break;
	  }
	  ++i;
	}
	proc_t result = i / nCells;
	return result;
      };
      hemelb::redblood::parallel::ExchangeCells xchange(graphComm);
      xchange.PostCellMessageLength(distributions, owned, ownership);
      xchange.PostCells(distributions, owned, ownership);
      auto const distCells = xchange.ReceiveCells(templates);
      auto const &lentCells = std::get<2>(distCells);

      // Actually perform velocity integration
      hemelb::redblood::parallel::IntegrateVelocities integrator(graphComm);
      integrator.PostMessageLength(lentCells);
      integrator.ComputeLocalVelocitiesAndUpdatePositions<Traits>(fieldData, owned);
      integrator.PostVelocities<Traits>(fieldData, lentCells);
      integrator.UpdatePositionsNonLocal(distributions, owned);

      if (world.Rank() == 0) {
	REQUIRE( (checkMoves - cells[0]->GetVertices()[0]).GetMagnitude() > 1e-8);
      }

      for (auto const &cell : cells) {
	using Positions = std::vector<LatticePosition>;
	Positions positions = world.Rank() == 0 ?
	  cell->GetVertices() :
	  Positions(cell->GetNumberOfNodes(), LatticePosition::Zero());
	world.Broadcast(to_span(positions), 0);
	if (world.Rank() != 0 and owned.find(cell) != owned.end()) {
	  for (auto const item : util::zip(cell->GetVertices(), positions)) {
	    REQUIRE(ApproxV(std::get<1>(item)).Margin(1e-8) == std::get<0>(item));
	  }
	}
      }
    }

    METHOD_AS_TEST_CASE(MPIIntegrateVelocitiesTests::testMidRegion<hemelb::redblood::stencil::FourPoint>,
			"Test velocity integration in mid region with 4 point stencil",
			"[redblood]");
    METHOD_AS_TEST_CASE(MPIIntegrateVelocitiesTests::testEdgeRegion<hemelb::redblood::stencil::FourPoint>,
			"Test velocity integration in edge region with 4 point stencil",
			"[redblood]");
    METHOD_AS_TEST_CASE(MPIIntegrateVelocitiesTests::testMidRegion<hemelb::redblood::stencil::FourPoint>,
			"Test velocity integration in all with 4 point stencil",
			"[redblood]");

    METHOD_AS_TEST_CASE(MPIIntegrateVelocitiesTests::testMidRegion<hemelb::redblood::stencil::ThreePoint>,
			"Test velocity integration in mid region with 3 point stencil",
			"[redblood]");
    METHOD_AS_TEST_CASE(MPIIntegrateVelocitiesTests::testEdgeRegion<hemelb::redblood::stencil::ThreePoint>,
			"Test velocity integration in edge region with 3 point stencil",
			"[redblood]");
    METHOD_AS_TEST_CASE(MPIIntegrateVelocitiesTests::testMidRegion<hemelb::redblood::stencil::ThreePoint>,
			"Test velocity integration in all with 3 point stencil",
			"[redblood]");

    METHOD_AS_TEST_CASE(MPIIntegrateVelocitiesTests::testMidRegion<hemelb::redblood::stencil::CosineApprox>,
			"Test velocity integration in mid region with cosine stencil",
			"[redblood]");
    METHOD_AS_TEST_CASE(MPIIntegrateVelocitiesTests::testEdgeRegion<hemelb::redblood::stencil::CosineApprox>,
			"Test velocity integration in edge region with cosine stencil",
			"[redblood]");
    METHOD_AS_TEST_CASE(MPIIntegrateVelocitiesTests::testMidRegion<hemelb::redblood::stencil::CosineApprox>,
			"Test velocity integration in all with cosine stencil",
			"[redblood]");

    METHOD_AS_TEST_CASE(MPIIntegrateVelocitiesTests::testMidRegion<hemelb::redblood::stencil::TwoPoint>,
			"Test velocity integration in mid region with 2 point stencil",
			"[redblood]");
    METHOD_AS_TEST_CASE(MPIIntegrateVelocitiesTests::testEdgeRegion<hemelb::redblood::stencil::TwoPoint>,
			"Test velocity integration in edge region with 2 point stencil",
			"[redblood]");
    METHOD_AS_TEST_CASE(MPIIntegrateVelocitiesTests::testMidRegion<hemelb::redblood::stencil::TwoPoint>,
			"Test velocity integration in all with 2 point stencil",
			"[redblood]");
}
