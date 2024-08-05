// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <map>
#include <catch2/catch.hpp>
#include "tests/helpers/FolderTestFixture.h"

#include "geometry/LookupTree.h"
#include "geometry/GeometryReader.h"
#include "lb/lattices/D3Q15.h"

namespace hemelb::tests
{
    using namespace geometry::octree;
    TEST_CASE("LookupTree - 1D bit twiddling", "[geometry]") {
        // Need to take 16 bit unsigned, spread the bits out into triples,
        // do anything to the "extra" bits and then go back ok.
        //
        // I.e.: abcd efgh ijkl mnop ->
        // 0000 0000 0000 0000 00a0 0b00 c00d 00e0 0f00 g00h 00i0 0j00 k00l 00m0 0n00 o00p
        std::uint16_t const i = GENERATE(0, 1, 2, 4, 8, 0xf, 0xaaaa, 0xffff);
        auto j = bitspread_one(i);
        static_assert(std::is_same_v<std::uint64_t, decltype(j)>);

        constexpr std::uint64_t HIGHMASK = 0xffff000000000000;
        REQUIRE((j & HIGHMASK) == 0);

        constexpr std::uint64_t EMPTYMASK = 0x0000db6db6db6db6;
        REQUIRE((j & EMPTYMASK) == 0);

        auto k = bitcondense_one(j ^ EMPTYMASK);
        static_assert(std::is_same_v<std::uint16_t, decltype(k)>);
        REQUIRE(k == i);
    }

    TEST_CASE("LookupTree - 3D bit twiddling", "[geometry]") {
        // Now interleave three 1D twiddled values and get them back
        std::uint16_t const i = GENERATE(0, 1, 0xf, 0xffff);
        std::uint16_t const j = GENERATE(0, 2, 0xf0, 0xffff);
        std::uint16_t const k = GENERATE(0, 4, 0xf00, 0xffff);

        auto ijk = Vec16(i,j,k);
        auto oct = ijk_to_oct(ijk);

        if (ijk == Vec16::Zero())
            REQUIRE(oct == 0U);
        else if (ijk == Vec16(0xffff))
            REQUIRE(oct == 0x0000ffffffffffffUL);

        auto xyz = oct_to_ijk(oct);
        REQUIRE(xyz[0] == i);
        REQUIRE(xyz[1] == j);
        REQUIRE(xyz[2] == k);
    }

    // Let's think about an example geometry that fits in a
    // 4-level tree (16 along each side) and is a simple cube.
    // These are the coords of the corners (usual up-to-but-not-including convention).
    struct SimpleTreeFixture {
        static constexpr auto const lo = Vec16(1, 3, 5);
        static constexpr auto hi = Vec16(9, 10, 11);
        static constexpr site_t block_size = 8;
        static constexpr site_t sites_per_block = 512;
        static constexpr site_t n_fluid_blocks = 8 * 7 * 6;

        LookupTree ref_tree;

        SimpleTreeFixture() : ref_tree(4) {
            // Build the tree, sort of by hand.
            for (auto block_ijk: IterBounds{hi}) {
                if (block_ijk[0] < lo[0] || block_ijk[1] < lo[1] || block_ijk[2] < lo[2])
                    continue;
                auto oct = ijk_to_oct(block_ijk);
                // leaf
                ref_tree.levels[4].node_ids.push_back(oct);
                ref_tree.levels[4].sites_per_node.push_back(sites_per_block);
                // next
                auto child_local = oct & 7U;
                auto parent_oct = oct >> 3;
                if (ref_tree.levels[3].node_ids.empty() || ref_tree.levels[3].node_ids.back() != parent_oct) {
                    ref_tree.levels[3].node_ids.push_back(parent_oct);
                    ref_tree.levels[3].sites_per_node.push_back(0);
                    ref_tree.levels[3].child_indices.push_back(Level::NOCHILDREN);
                }
                ref_tree.levels[3].sites_per_node.back() += sites_per_block;
                auto &link_34 = ref_tree.levels[3].child_indices.back()[child_local];
                if (link_34 == Level::NC) {
                    link_34 = ref_tree.levels[4].node_ids.size() - 1;
                } else {
                    REQUIRE(link_34 == (ref_tree.levels[4].node_ids.size() - 1));
                }

                oct = parent_oct;
                child_local = oct & 7U;
                parent_oct = oct >> 3;
                if (ref_tree.levels[2].node_ids.empty() || ref_tree.levels[2].node_ids.back() != parent_oct) {
                    ref_tree.levels[2].node_ids.push_back(parent_oct);
                    ref_tree.levels[2].sites_per_node.push_back(0);
                    ref_tree.levels[2].child_indices.push_back(Level::NOCHILDREN);
                }
                ref_tree.levels[2].sites_per_node.back() += sites_per_block;
                auto &link_23 = ref_tree.levels[2].child_indices.back()[child_local];
                if (link_23 == Level::NC) {
                    link_23 = ref_tree.levels[3].node_ids.size() - 1;
                } else {
                    REQUIRE(link_23 == (ref_tree.levels[3].node_ids.size() - 1));
                }

                oct = parent_oct;
                child_local = oct & 7U;
                parent_oct = oct >> 3;
                if (ref_tree.levels[1].node_ids.empty() || ref_tree.levels[1].node_ids.back() != parent_oct) {
                    ref_tree.levels[1].node_ids.push_back(parent_oct);
                    ref_tree.levels[1].sites_per_node.push_back(0);
                    ref_tree.levels[1].child_indices.push_back(Level::NOCHILDREN);
                }
                ref_tree.levels[1].sites_per_node.back() += sites_per_block;
                auto &link_12 = ref_tree.levels[1].child_indices.back()[child_local];
                if (link_12 == Level::NC) {
                    link_12 = ref_tree.levels[2].node_ids.size() - 1;
                } else {
                    REQUIRE(link_12 == (ref_tree.levels[2].node_ids.size() - 1));
                }

                oct = parent_oct;
                child_local = oct & 7U;
                parent_oct = oct >> 3;
                REQUIRE(parent_oct == 0);
                if (ref_tree.levels[0].node_ids.empty() || ref_tree.levels[0].node_ids.back() != parent_oct) {
                    ref_tree.levels[0].node_ids.push_back(parent_oct);
                    ref_tree.levels[0].sites_per_node.push_back(0);
                    ref_tree.levels[0].child_indices.push_back(Level::NOCHILDREN);
                }
                ref_tree.levels[0].sites_per_node.back() += sites_per_block;
                auto &link_01 = ref_tree.levels[0].child_indices.back()[child_local];
                if (link_01 == Level::NC) {
                    link_01 = ref_tree.levels[1].node_ids.size() - 1;
                } else {
                    REQUIRE(link_01 == (ref_tree.levels[1].node_ids.size() - 1));
                }

            }
        }
    };

    TEST_CASE_METHOD(SimpleTreeFixture, "LookupTree - test iterator", "[geometry]") {
        std::vector<Vec16> blocks;
        std::vector<std::uint64_t> oct_ids;
        std::uint64_t last_oct;
        for (auto block_ijk: IterBounds{hi}) {
            auto oct = ijk_to_oct(block_ijk);
            // octree order
            if (!blocks.empty())
                REQUIRE(oct > last_oct);

            blocks.push_back(block_ijk);
            last_oct = oct;

            // In range
            for (auto d: {0,1,2}) {
                // This test is not correct as want to consider every potential block in the domain.
                // REQUIRE(block_ijk[d] >= lo[d]);
                REQUIRE(block_ijk[d] < hi[d]);
            }
        }
        REQUIRE(blocks.size() == hi[0]*hi[1]*hi[2]);
    }

    TEST_CASE_METHOD(SimpleTreeFixture, "LookupTree - check SimpleTreeFixture", "[geometry]") {
        REQUIRE(ref_tree.levels[0].node_ids.size() == 1);
        REQUIRE(ref_tree.levels[0].sites_per_node.size() == 1);
        REQUIRE(ref_tree.levels[0].sites_per_node[0] == n_fluid_blocks * sites_per_block);
        REQUIRE(ref_tree.levels[4].node_ids.size() == n_fluid_blocks);
        REQUIRE(ref_tree.levels[4].sites_per_node.size() == n_fluid_blocks);
    }

    template <typename T>
    void check_vec(std::vector<T> const& actual, std::vector<T> const& expected) {
        REQUIRE(actual.size() == expected.size());
        for (int i = 0; i < std::ssize(actual); ++i) {
            REQUIRE(actual[i] == expected[i]);
        }
    }
    void check_level(Level const& actual, Level const& expected) {
        REQUIRE(actual.level == expected.level);
        check_vec(actual.node_ids, expected.node_ids);
        check_vec(actual.sites_per_node, expected.sites_per_node);
        check_vec(actual.child_indices, actual.child_indices);
    }

    TEST_CASE_METHOD(SimpleTreeFixture, "LookupTree - builder works for simple case", "[geometry]") {
        std::vector<site_t> fluidSitesPerBlock(990);
        int ijk = 0;
        for (int i = 0; i < hi[0]; ++i)
            for (int j = 0; j < hi[1]; ++j)
                for (int k = 0; k < hi[2]; ++k) {
                    fluidSitesPerBlock[ijk] = (i < lo[0] || j < lo[1] || k < lo[2]) ? 0 : sites_per_block;
                    ++ijk;
                }
        REQUIRE(ijk == 990);

        LookupTree actual = build_block_tree(hi, fluidSitesPerBlock);

        REQUIRE(actual.n_levels == ref_tree.n_levels);

        for (U16 i = 0; i < actual.n_levels; ++i) {
            check_level(actual.levels[i], ref_tree.levels[i]);
        }
    }

    TEST_CASE_METHOD(helpers::FolderTestFixture, "LookupTree - build tree from GMY", "[geometry]") {
        // CopyResourceToTempdir("large_cylinder.xml");
        CopyResourceToTempdir("large_cylinder.gmy");
        MoveToTempdir();

        auto timings = std::make_unique<reporting::Timers>();
        auto reader = std::make_unique<geometry::GeometryReader>(lb::D3Q15::GetLatticeInfo(),
                                                                 *timings,
                                                                 Comms());
        auto result = reader->LoadAndDecompose("large_cylinder.gmy");
        auto& tree = result.block_store->GetTree();
        //auto tree = //result.StealBlockTree();
        // hlb-gmy-countsites -v -v gives:
        //  BlockCounts = [2 2 5]
        //  BlockSize = 8
        //  TotalFluidSites = 5576
        //  BlocksWithFluidSites = 20
        //    Index  Sites Oct ID
        //    0  287  0
        //    1  328  1
        //    2  328  8
        //    3  328  9
        //    4  123  64
        //    5  287  2
        //    6  328  3
        //    7  328  10
        //    8  328  11
        //    9  123  66
        //    10  287  4
        //    11  328  5
        //    12  328  12
        //    13  328  13
        //    14  123  68
        //    15  287  6
        //    16  328  7
        //    17  328  14
        //    18  328  15
        //    19  123  70
        REQUIRE(tree.n_levels == 3);
        REQUIRE(tree.levels[0].sites_per_node.size() == 1);
        REQUIRE(tree.levels[0].sites_per_node[0] == 5576);
        REQUIRE(tree.levels[3].node_ids.size() == 20);

        auto b000 = tree.GetPath({0,0,0});
        REQUIRE(b000.path[0] == 0);

        auto b114 = tree.GetPath({1,1,4});
        // 1,1,4 => 0001, 0001, 0100 => 000 001 000 110
        REQUIRE(tree.levels[0].node_ids[b114.path[0]] == 0);
        REQUIRE(tree.levels[1].node_ids[b114.path[1]] == 1);
        REQUIRE(tree.levels[2].node_ids[b114.path[2]] == 8);
        REQUIRE(tree.levels[3].node_ids[b114.path[3]] == 70);

        auto b334 = tree.GetPath({3,3,4});
        // 3,3,4 => 0011, 0011, 0100 => 000 001 110 110
        REQUIRE(tree.levels[0].node_ids[b334.path[0]] == 0);
        REQUIRE(tree.levels[1].node_ids[b334.path[1]] == 1);
        REQUIRE(b334.path[2] == Level::NC);
        REQUIRE(b334.path[3] == Level::NC);

        // fluid sites per block, ordered by GMY index
        auto spb = std::vector<U64>{
            287, 328, 328, 328, 123, 287, 328, 328, 328, 123, 287, 328, 328, 328, 123, 287, 328, 328, 328, 123
        };
        std::map<U64, U64> oct2spb;
        for (U64 block = 0U; block < 20; ++block) {
            auto b_ijk = Vec16(block / 10U, (block / 5U) % 2U, block % 5U);
            auto oct = ijk_to_oct(b_ijk);
            oct2spb[oct] = spb[block];
        }

        for (int i = 0; i < 20; ++i) {
            auto oct = tree.levels[3].node_ids[i];
            REQUIRE(oct2spb.contains(oct));
            REQUIRE(tree.levels[3].sites_per_node[i] == oct2spb[oct]);
        }

        // The next level up in the tree, have three nodes
        REQUIRE(tree.levels[2].node_ids.size() == 3);
        REQUIRE(tree.levels[2].node_ids[0] == 0);
        REQUIRE(tree.levels[2].node_ids[1] == 1);
        REQUIRE(tree.levels[2].node_ids[2] == 8);
        REQUIRE(tree.levels[2].sites_per_node[0] == 2460);
        REQUIRE(tree.levels[2].sites_per_node[1] == 2624);
        REQUIRE(tree.levels[2].sites_per_node[2] == 492);

        // The next level up in the tree, have two nodes
        REQUIRE(tree.levels[1].node_ids.size() == 2);
        REQUIRE(tree.levels[1].node_ids[0] == 0);
        REQUIRE(tree.levels[1].node_ids[1] == 1);
        REQUIRE(tree.levels[1].sites_per_node[0] == 5084);
        REQUIRE(tree.levels[1].sites_per_node[1] == 492);
    }
}
