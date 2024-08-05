// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <iterator>
#include <catch2/catch.hpp>
#include "redblood/WallCellPairIterator.h"
#include "lb/lattices/D3Q15.h"

#include "tests/helpers/FourCubeLatticeData.h"
#include "tests/helpers/HasCommsTestFixture.h"
#include "tests/helpers/SiteIterator.h"
#include "tests/helpers/LatticeDataAccess.h"

namespace hemelb::tests {
    using namespace redblood;

    TEST_CASE_METHOD(helpers::HasCommsTestFixture, "WallCellPairIteratorTests", "[redblood]") {
        LatticeDistance const cutoff = 3.0;
        LatticeDistance const interactionDistance = 0.5;
        LatticeDistance const halo = interactionDistance + 1e-6;
        using Lattice = lb::D3Q15;
        std::unique_ptr<tests::FourCubeLatticeData> latticeData{FourCubeLatticeData::Create(Comms(), 27 + 2)};
        auto &dom = latticeData->GetDomain();

        for (auto const site: std::cref(*latticeData)) {
            if (!site.IsWall())
                continue;

            for (Direction d(0); d < Lattice::NUMVECTORS; ++d) {
                if (site.HasWall(d)) {
                    dom.SetBoundaryDistance(site.GetIndex(), d, 0.5);
                }
            }
        }

        SECTION("testOneCellNode") {
            auto const wallDnC = createWallNodeDnC<Lattice>(dom, cutoff, halo);
            auto const cell = std::make_shared<Cell>(tetrahedron());
            *cell *= 3;
            *cell += LatticePosition{100};
            cell->GetVertices()[0] = LatticePosition(0.5, 3.5 * cutoff, 3.5 * cutoff);

            DivideConquerCells const cellDnC({cell}, cutoff, interactionDistance);

            auto range = iterate(cellDnC, wallDnC, interactionDistance);
            auto const last_iterator = std::end(range);
            auto const first_iterator = std::begin(range);
            REQUIRE(first_iterator != last_iterator);
            REQUIRE(std::ptrdiff_t(4) == std::distance(first_iterator, last_iterator));
            for (auto const item: range) {
                REQUIRE((item.wallNode - item.cellNode).GetMagnitude() < interactionDistance);
            }
            auto short_range = iterate(cellDnC, wallDnC, 0.05);
            REQUIRE(std::ptrdiff_t(0) == std::distance(
                    std::begin(short_range),
                    std::end(short_range)
            ));
        }

        SECTION("testNoWall") {
            DivideConquer<WallNode> const wallDnC(3e0);
            auto const cell = std::make_shared<Cell>(tetrahedron());
            *cell *= 3;
            *cell += LatticePosition{100};
            cell->GetVertices()[0] = LatticePosition(0.5, 3.5 * cutoff, 3.5 * cutoff);
            DivideConquerCells const cellDnC({cell}, cutoff, interactionDistance);

            auto range = iterate(cellDnC, wallDnC, interactionDistance);
            auto const last_iterator = std::end(range);
            auto const first_iterator = std::begin(range);
            REQUIRE(first_iterator == last_iterator);
            REQUIRE(std::ptrdiff_t(0) == std::distance(first_iterator, last_iterator));
        }

        SECTION("testTwoCellNodes") {
            auto const wallDnC = createWallNodeDnC<Lattice>(dom, cutoff, halo);
            auto const cell = std::make_shared<Cell>(tetrahedron());
            *cell *= 3;
            *cell += LatticePosition{100};
            cell->GetVertices()[0] = LatticePosition(0.5, 3.5 * cutoff, 3.5 * cutoff);
            cell->GetVertices()[1] = LatticePosition(0.49, 3.5 * cutoff, 3.5 * cutoff);

            DivideConquerCells const cellDnC({cell}, cutoff, interactionDistance);

            auto range = iterate(cellDnC, wallDnC, interactionDistance);
            REQUIRE(std::ptrdiff_t(8) == std::distance(std::begin(range), std::end(range)));
        }

        SECTION("testHaloAndNeighboringBoxes") {
            auto testHaloAndNeighboringBoxes = [&](Dimensionless where, std::ptrdiff_t howMany) {
                auto const wallDnC = createWallNodeDnC<Lattice>(dom, cutoff, halo);
                auto const cell = std::make_shared<Cell>(tetrahedron());
                *cell *= 3;
                *cell += LatticePosition{100};
                cell->GetVertices()[0] = LatticePosition(0.5, where * cutoff, where * cutoff);

                DivideConquerCells const cellDnC({cell}, cutoff, interactionDistance);
                auto range = iterate(cellDnC,
                                     wallDnC,
                                     interactionDistance);
                REQUIRE(howMany == std::distance(std::begin(range), std::end(range)));
            };

            // center of a divide and conquer box and in the middle of a square of fluid-sites
            testHaloAndNeighboringBoxes(3.5, 4);
            // center of a divide and conquer box and on top of a fluid-sites
            testHaloAndNeighboringBoxes(3.33333333, 5);
            // // close to range border and in the middle of a square of fluid-sites
            testHaloAndNeighboringBoxes(3.0, 5);
        }

        auto testAllWallNodesFound = [&](Dimensionless where) {
            auto const wallDnC = createWallNodeDnC<Lattice>(dom, cutoff, halo);
            auto const cell = std::make_shared<Cell>(tetrahedron());
            *cell *= 3;
            *cell += LatticePosition{100};
            cell->GetVertices()[0] = LatticePosition(0.5, where * cutoff, where * cutoff);
            DivideConquerCells const cellDnC({cell}, cutoff, interactionDistance);

            // Creates list of wall-nodes and check unicity
            constexpr double TOL = 1e-8;

            std::vector<LatticePosition> wallNodes;
            for (auto const nodes: iterate(cellDnC, wallDnC, interactionDistance)) {
                auto const same = std::find_if(wallNodes.begin(),
                                               wallNodes.end(),
                                               [&](LatticePosition const &x) {
                                                   return (nodes.wallNode - x).GetMagnitude() < TOL;
                                               });
                REQUIRE(same == wallNodes.end());
                wallNodes.push_back(nodes.wallNode);
            }

            // Loop over all nodes, check whether they are in range, check they are in list
            for (auto const site: std::cref(*latticeData)) {
                if (!site.IsWall())
                    continue;

                for (Direction d = 0; d < Lattice::NUMVECTORS; ++d) {
                    if (!site.HasWall(d))
                        continue;
                    auto const node =
                            Lattice::CD[d].GetNormalised() * site.GetWallDistance<Lattice>(d)
                            + site.GetGlobalSiteCoords();
                    auto const visited = std::find_if(wallNodes.begin(),
                                                      wallNodes.end(),
                                                      [&](LatticePosition const &x) {
                                                          return (node - x).GetMagnitude() < TOL;
                                                      });
                    bool inRange = (node - cell->GetVertices()[0]).GetMagnitude() < interactionDistance;
                    REQUIRE((visited != wallNodes.end()) == inRange);
                    if (visited != wallNodes.end()) {
                        wallNodes.erase(visited);
                    }
                }
            }
            // At this point both methods found the same nodes if and only if wallNodes is empty
            REQUIRE(size_t(0) == wallNodes.size());
        };

        SECTION("testAllWallNodesFound") {
            size_t const N = 10;
            for (size_t i(0); i <= N; ++i) {
                testAllWallNodesFound(Dimensionless(i) / Dimensionless(2 * N) + 3.0);
                testAllWallNodesFound(Dimensionless(i) / Dimensionless(2 * N) + 0.0);
            }
        }

        SECTION("testInteractionPassedOnToFluid") {
            using STENCIL = stencil::FourPoint;
            size_t const N = 10;
            for (size_t i(0); i <= N; ++i) {
                //testInteractionPassedOnToFluid<STENCIL>();
                Dimensionless where = Dimensionless(i) / Dimensionless(2 * N) + 3.0;
                auto const wallDnC = createWallNodeDnC<Lattice>(dom, cutoff, halo);
                auto const cell = std::make_shared<Cell>(tetrahedron());
                LatticePosition const node(0.6, where * cutoff, where * cutoff);
                *cell *= 3;
                *cell += LatticePosition{100};
                cell->GetVertices()[0] = node;
                DivideConquerCells const cellDnC({cell}, cutoff, interactionDistance);

                // Set forces to zero
                helpers::ZeroOutForces(latticeData.get());

                // Finds pairs, computes interaction, spread forces to lattice
                addCell2WallInteractions<STENCIL>(DivideConquerCells({cell}, cutoff, halo),
                                                  wallDnC,
                                                  Node2NodeForce(1.0, interactionDistance),
                                                  *latticeData);

                for (auto const site: std::cref(*latticeData)) {
                    auto const d = LatticePosition(site.GetGlobalSiteCoords()) - node;
                    REQUIRE((STENCIL::stencil(d) > 1e-12) == (site.GetForce().GetMagnitude() > 1e-12));
                }

            }
        }
    }
}
