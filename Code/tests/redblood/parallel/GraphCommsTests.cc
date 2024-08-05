// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "tests/redblood/Fixtures.h"
#include "tests/redblood/parallel/ParallelFixtures.h"
#include "lb/iolets/InOutLet.h"

namespace hemelb::tests
{
    using namespace redblood;

    class GraphCommsTests : public OpenSimFixture
    {
    public:
        GraphCommsTests();

        void testOrderingOfVectors();
        void testDumbGraphCommunicator();
        void testGraphCommunicator();
        void testComputeCellsEffectiveSize();
        void testComputeGlobalCoordsToProcMap();
    };

    GraphCommsTests::GraphCommsTests() : OpenSimFixture() {
      // Have everything ready to creates simulations
      if (Comms().Rank() == 0) {
	CopyResourceToTempdir("large_cylinder_rbc.xml");
	CopyResourceToTempdir("large_cylinder.gmy");
	CopyResourceToTempdir("red_blood_cell.txt");

	// This simulation duration is sufficient to pick up the original force spreading
	// issue that motivated the test. Run the test for longer in order to check other
	// aspects of the parallel implementation against a sequential run.
	ModifyXMLInput("large_cylinder_rbc.xml", { "simulation", "steps", "value" }, 1000);
      }
      Comms().Barrier();

      options =
	std::make_shared<configuration::CommandLine>(configuration::CommandLine { "hemelb",
	      "-in",
	      "large_cylinder_rbc.xml",
    });
    }

    void GraphCommsTests::testOrderingOfVectors() {
      auto cmp = parallel::VectorOctreeOrdering{};

      Vec16 foo{0, 0, 0};
      Vec16 bar{1, 0, 0};

      // // Catch implementation
      // CHECK(cmp(foo, bar));

      // bar = {0, 1, 0};
      // CHECK(cmp(foo, bar));

      // bar = {0, 0, 1};
      // CHECK(cmp(foo, bar));

      // bar = {0, 0, 0};
      // CHECK(!cmp(foo, bar));

      // foo = {0, 0, 1};
      // CHECK(!cmp(foo, bar));
      REQUIRE(cmp(foo, bar));

      bar = {0, 1, 0};
      REQUIRE(cmp(foo, bar));

      bar = {0, 0, 1};
      REQUIRE(cmp(foo, bar));

      bar = {0, 0, 0};
      REQUIRE(!cmp(foo, bar));

      foo = {0, 0, 1};
      REQUIRE(!cmp(foo, bar));
    }

    void GraphCommsTests::testDumbGraphCommunicator()
    {
        auto graph = CreateDumbGraphComm(Comms());
        auto my_neighbours = graph.GetNeighbors();

        REQUIRE(std::ssize(my_neighbours) == Comms().Size() - 1);

        // All other processes must be neighbours
        for (int neighbour = 0; neighbour < Comms().Size(); ++neighbour) {
            if (neighbour != Comms().Rank()) {
                REQUIRE(
                        std::find(my_neighbours.begin(), my_neighbours.end(), neighbour) != my_neighbours.end()
                );
            }
        }
    }

    void GraphCommsTests::testGraphCommunicator()
    {
        using parallel::ComputeProcessorNeighbourhood;

        // Test only makes sense if run with 4 cores
        auto comms = Comms();
        REQUIRE(4 == comms.Size());

        // Setup simulation with cylinder
        auto sim = CreateSim<stencil::FourPoint>(comms);
        REQUIRE(sim);
        auto &latticeData = sim->GetFieldData();

        // Compute neighbourhoods (cylinder is 4.8e-5 long, a cell
        // effective size of 2e-6 won't let cells span across more than
        // two subdomains)
        auto neighbourhoods =
                ComputeProcessorNeighbourhood(comms,
                                              latticeData.GetDomain(),
                                              2e-6 / sim->GetSimConfig().GetVoxelSize());

        // Parmetis divides the cylinder in four consecutive cylindrical
        // subdomains with interfaces roughly parallel to the iolets.
        // The ranks are ordered 0,1,2,3 along the positive direction of
        // the z axis.
        std::vector<std::vector<int>> expected_neighbourhoods;
        expected_neighbourhoods.push_back( { 1 });
        expected_neighbourhoods.push_back( { 0, 2 });
        expected_neighbourhoods.push_back( { 1, 3 });
        expected_neighbourhoods.push_back( { 2 });

        // Compare computed vs expected neighbourhood for this rank
        REQUIRE(expected_neighbourhoods[comms.Rank()] == neighbourhoods);
    }

    void GraphCommsTests::testComputeCellsEffectiveSize() {
        using parallel::ComputeCellsEffectiveSize;

        // Setup simulation with cylinder
        auto comms = Comms();
        auto sim = CreateSim<stencil::FourPoint>(comms);
        REQUIRE(sim);

        auto simConf = sim->GetSimConfig();
        REQUIRE(simConf.HasRBCSection());
        auto builder = configuration::SimBuilder(simConf);
        auto rbcConf = simConf.GetRBCConfig();

        auto meshes = [&]() {
            using IoletPtr = util::clone_ptr<lb::InOutLet>;
            auto ccb = redblood::CellControllerBuilder(builder.GetUnitConverter());
            auto inlets = builder.BuildIolets(simConf.GetInlets());
            auto outlets = builder.BuildIolets(simConf.GetOutlets());
            auto mk_view = [](std::vector<IoletPtr> const &iolets) {
                return redblood::CountedIoletView(
                        [&iolets]() { return iolets.size(); },
                        [&iolets](unsigned i) { return iolets[i].get(); }
                );
            };
            return ccb.build_template_cells(rbcConf.meshes, mk_view(inlets), mk_view(outlets));
        }();

        // Biggest cell radius in lattice units times a tolerance
        REQUIRE(
                Approx(
                        parallel::MAXIMUM_SIZE_TO_RADIUS_RATIO * (8e-06 / simConf.GetVoxelSize())
                ).margin(1e-9) == ComputeCellsEffectiveSize(*meshes)
        );
    }

    void GraphCommsTests::testComputeGlobalCoordsToProcMap()
    {
        using parallel::ComputeGlobalCoordsToProcMap;
        using parallel::ComputeProcessorNeighbourhood;

        // Test only makes sense if run with 4 cores
        auto comms = Comms();
        REQUIRE(4 == comms.Size());

        // Setup simulation with cylinder
        auto sim = CreateSim<stencil::FourPoint>(comms);
        REQUIRE(sim);
        auto &domain = sim->GetFieldData().GetDomain();
        auto graphComm = comms.DistGraphAdjacent(
              ComputeProcessorNeighbourhood(comms,
                                            domain,
                                            2e-6 / sim->GetSimConfig().GetVoxelSize())
        );
        auto const& globalCoordsToProcMap = ComputeGlobalCoordsToProcMap(graphComm, domain);

        // The first lattice site for each rank is
        std::vector<Vec16> lattice_coords;
        lattice_coords.emplace_back(2, 6, 2);
        lattice_coords.emplace_back(2, 6, 9);
        lattice_coords.emplace_back(2, 6, 19);
        lattice_coords.emplace_back(2, 6, 29);

        // Parmetis divides the cylinder in four consecutive cylindrical subdomains.
        // The ranks are ordered 0,1,2,3 along the positive direction of the z axis.
        // The cell size is such that only ranks to the left and right are neighbours.
        std::vector<std::vector<proc_t>> neighs_to_know_about;
        neighs_to_know_about.push_back({0, 1});
        neighs_to_know_about.push_back({0, 1, 2});
        neighs_to_know_about.push_back({1, 2, 3});
        neighs_to_know_about.push_back({2, 3});

        // Each process should know who is the owner of a local lattice and a lattice belonging to each neighbour
        for (auto neigh : neighs_to_know_about[comms.Rank()]) {
            // Use map::at to throw when the lattice coordinates are not known
            REQUIRE(globalCoordsToProcMap.at(lattice_coords[neigh]) == neigh);
        }
    }

    METHOD_AS_TEST_CASE(GraphCommsTests::testOrderingOfVectors,
			"testOrderingOfVectors",
			"[redblood]");
    METHOD_AS_TEST_CASE(GraphCommsTests::testDumbGraphCommunicator,
			"testDumbGraphCommunicator",
			"[redblood]");
    METHOD_AS_TEST_CASE(GraphCommsTests::testGraphCommunicator,
			"testGraphCommunicator",
			"[redblood]");
    METHOD_AS_TEST_CASE(GraphCommsTests::testComputeCellsEffectiveSize,
			"testComputeCellsEffectiveSize",
			"[redblood]");
    METHOD_AS_TEST_CASE(GraphCommsTests::testComputeGlobalCoordsToProcMap,
			"testComputeGlobalCoordsToProcMap",
			"[redblood]");
}
