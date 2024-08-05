// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <cstdio>
#include <unistd.h>

#include <catch2/catch.hpp>

#include "Traits.h"
#include "configuration/SimBuilder.h"
#include "redblood/Cell.h"
#include "redblood/CellController.h"
#include "tests/redblood/Fixtures.h"

#include "tests/helpers/ApproxVector.h"
#include "tests/helpers/LatticeDataAccess.h"
#include "tests/helpers/FolderTestFixture.h"

namespace hemelb::tests
{
    using namespace redblood;

    template <typename STENCIL>
    class NodeIntegrationTestsFixture : public helpers::FolderTestFixture {
        using TRAITS = Traits<
                lb::DefaultLattice, lb::GuoForcingLBGK, lb::Normal,
                lb::DefaultStreamer, lb::DefaultWallStreamer, lb::DefaultInletStreamer, lb::DefaultOutletStreamer,
                STENCIL
        >;
    protected:
        static constexpr char const* xml_name = "large_cylinder_rbc.xml";
        static constexpr int argc = 3;
        static const char* argv[argc];

    public:
        NodeIntegrationTestsFixture() {
            this->CopyResourceToTempdir("red_blood_cell.txt");
            CopyResourceToTempdir("large_cylinder.gmy");
        }

        //! Creates a simulation. Does not run it.
        auto createSim(size_t steps,
                                                                   Dimensionless cell,
                                                                   Dimensionless wall) const {
            CopyResourceToTempdir(xml_name);
            if (std::filesystem::exists("results"))
            {
                std::filesystem::remove_all("results");
            }
            DeleteXMLInput(xml_name, { "inlets", "inlet", "insertcell" });
            DeleteXMLInput(xml_name, { "inlets", "inlet", "flowextension" });
            DeleteXMLInput(xml_name, { "outlets", "outlet", "flowextension" });
            DeleteXMLInput(xml_name, { "properties", "propertyoutput" });
            ModifyXMLInput(xml_name, { "inlets", "inlet", "condition", "mean", "value" }, 0);
            ModifyXMLInput(xml_name, { "simulation", "steps", "value" }, steps);
            ModifyXMLInput(xml_name, { "redbloodcells", "cell2Cell", "intensity", "value" }, cell);
            ModifyXMLInput(xml_name, { "redbloodcells", "cell2Cell", "cutoffdistance", "value" }, 1);
            ModifyXMLInput(xml_name, { "redbloodcells", "cell2Wall", "intensity", "value" }, wall);
            ModifyXMLInput(xml_name, { "redbloodcells", "cell2Wall", "cutoffdistance", "value" }, 1);
            auto options = std::make_shared<configuration::CommandLine>(argc, argv);
            auto sim = configuration::SimBuilder::CreateSim<TRAITS>(*options, Comms());
            helpers::LatticeDataAccess(&sim->GetFieldData()).ZeroOutForces();
            return sim;
        }

        //! Runs simulation with a single node
        LatticePosition testNodeWall(Dimensionless intensity, LatticePosition const &where) {

            using CellController = CellController<TRAITS>;

            auto const sim = this->createSim(1, 0, intensity);
            auto controller = std::static_pointer_cast<CellController>(sim->GetCellController());
            REQUIRE(controller); // XML file contains RBC problem definition, therefore a CellController should exist already
            controller->AddCell(std::make_shared<NodeCell>(where));

            sim->RunSimulation();
            return (*controller->GetCells().begin())->GetVertices().front();
        }

        std::pair<LatticePosition, LatticePosition>
        testNodeNode(
                Dimensionless intensity, LatticePosition const & n0, LatticePosition const & n1
        ) {
            using CellController = hemelb::redblood::CellController<TRAITS>;

            auto const sim = this->createSim(3, intensity, 0);
            auto controller = std::static_pointer_cast<CellController>(sim->GetCellController());
            REQUIRE(controller); // XML file contains RBC problem definition, therefore a CellController should exist already
            auto const firstCell = std::make_shared<NodeCell>(n0);
            auto const secondCell = std::make_shared<NodeCell>(n1);
            controller->AddCell(firstCell);
            controller->AddCell(secondCell);

            sim->RunSimulation();
            return {
                    firstCell->GetVertices().front(),
                    secondCell->GetVertices().front()
            };
        }
    };
    template <typename STENCIL>
    const char* NodeIntegrationTestsFixture<STENCIL>::argv[NodeIntegrationTestsFixture<STENCIL>::argc] = {
            "hemelb", "-in", xml_name
    };

    TEMPLATE_TEST_CASE_METHOD(NodeIntegrationTestsFixture,
                              "A single node cell will repel from the wall",
                              "[redblood]",
                              stencil::FourPoint, stencil::CosineApprox, stencil::ThreePoint, stencil::TwoPoint) {
        auto approx = Approx(0.0).margin(1e-8);

        SECTION("A cell beyond the cutoff will not move") {
            auto const farFromWall = this->testNodeWall(2e0, { 7e0, 7e0, 16.1e0 });
            REQUIRE(approx(7e0) == farFromWall.x());
            REQUIRE(approx(7e0) == farFromWall.y());
            REQUIRE(approx(16.1e0) == farFromWall.z());
        }

        LatticePosition const nearWall(3.5, 2.5, 16.0);
        SECTION("no interaction because intensity is zero, even if near wall") {
            auto const noForce = this->testNodeWall(0e0, nearWall);
            REQUIRE(approx(nearWall.x()) == noForce.x());
            REQUIRE(approx(nearWall.y()) == noForce.y());
            REQUIRE(approx(nearWall.z()) == noForce.z());
        }

        SECTION("A cell near the wall will repel from it in a normal direction (if symmetrically placed)") {
            // yes we can
            auto const withForce = this->testNodeWall(1e0, nearWall);
            auto const dr = withForce - nearWall;
            REQUIRE(dr.x() > 1e-4);
            REQUIRE(dr.y() > 1e-4);
            REQUIRE(approx(dr.z()) == 0.0);
        }

        SECTION("As forces depend on node location this will not be symmetric in general") {
            auto const start = nearWall + LatticePosition(0, 0, 0.1);
            auto const nonSym = this->testNodeWall(1e0, start);
            auto const dr = nonSym - start;
            REQUIRE(dr.x() > 1e-4);
            REQUIRE(dr.y() > 1e-4);
            REQUIRE(dr.z() > 1e-4);
        }
    }

    TEMPLATE_TEST_CASE_METHOD(NodeIntegrationTestsFixture,
                              "Two single node cells will repel eachother",
                              "[redblood]",
                              stencil::FourPoint, stencil::CosineApprox, stencil::ThreePoint, stencil::TwoPoint) {
        auto approx = Approx(0.0).margin(1e-8);

        LatticePosition const center{7, 7, 20};
        LatticePosition const z{0, 0, 1};

        SECTION("no interaction because far from each other") {
            auto half_separation = 2e0;
            auto const a = center - z * half_separation;
            auto const b = center + z * half_separation;
            auto const distant = this->testNodeNode(2e0, a, b);
            REQUIRE(distant.first == ApproxV(a));
            REQUIRE(distant.second == ApproxV(b));
        }

        SECTION("no interaction because intensity is zero") {
            auto half_separation = 0.2;
            auto const a = center - z * half_separation;
            auto const b = center + z * half_separation;
            auto const noForce = this->testNodeNode(0e0, a, b);
            REQUIRE(noForce.first == ApproxV(a));
            REQUIRE(noForce.second == ApproxV(b));
        }
        SECTION("cells nearby repel, symmetrically") {
            auto half_separation = 0.2;
            auto const a = center - z * half_separation;
            auto const b = center + z * half_separation;
            auto const withForce = this->testNodeNode(1e0, a, b);
            INFO(withForce.first);
            INFO(withForce.second);
            // a moves in negative z direction
            auto const da = withForce.first - a;
            REQUIRE(da.x() == approx(0));
            REQUIRE(da.y() == approx(0));
            REQUIRE(da.z() < -1e-4);
            // b moves in positive z direction
            auto const db = withForce.second - b;
            REQUIRE(db.x() == approx(0));
            REQUIRE(db.y() == approx(0));
            REQUIRE(db.z() > +1e-4);
            // Symmetric movement
            REQUIRE(da.z() + db.z() == approx(0));
        }
    }
}
