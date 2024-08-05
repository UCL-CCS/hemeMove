// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "geometry/decomposition/BasicDecomposition.h"
#include "geometry/GmyReadResult.h"
#include "geometry/LookupTree.h"
#include "log/Logger.h"
#include "net/mpi.h"

namespace hemelb::geometry::decomposition
{

    BasicDecomposition::BasicDecomposition(const GmyReadResult& geometry,
                                           int s) :
            geometry(geometry), comm_size(s)
    {
    }

    // Figure out how to divide a range between two groups of ranks of almost equal size
    auto share_range(std::vector<U64>::const_iterator begin, std::vector<U64>::const_iterator end, int N) {
        auto lo = *begin;
        auto hi = *end;
        auto delta = hi - lo;

        auto n_lo = N / 2;
        auto mid = float(lo) + float(delta) * float(n_lo) / float(N);
        auto middle = std::lower_bound(begin, end, mid);
        if ((*middle - mid) / float(delta) > 0.5f)
            --middle;

        return std::make_pair(n_lo, middle);
    }

    // Assign a group of processes of size N (i.e. ranks 0.. N-1) to the blocks with cumulative site counts.
    // Recursively split the range in half and assign depth-first.
    void assign_range(std::vector<U64>::const_iterator count_begin, std::vector<U64>::const_iterator count_end,
                      std::vector<int>::iterator rank_begin, std::vector<int>::iterator rank_end, int N);
    void assign_range(std::vector<U64>::const_iterator count_begin, std::vector<U64>::const_iterator count_end,
                      std::vector<int>::iterator rank_begin, std::vector<int>::iterator rank_end, int N) {
        if (N < 2)
            return;

        // Fairly split the communicator in half(ish)
        auto [n_lo, count_middle] = share_range(count_begin, count_end, N);
        auto dn = std::distance(count_begin, count_middle);
        auto rank_middle = rank_begin + dn;
        // Recursively assign the two halves
        assign_range(count_begin, count_middle, rank_begin, rank_middle, n_lo);
        assign_range(count_middle, count_end, rank_middle, rank_end, N - n_lo);
        // The second half needs to have the ranks incremented
        for (; rank_middle != rank_end; ++rank_middle)
            *rank_middle += n_lo;
    }

    std::vector<int> BasicDecomposition::Decompose(octree::LookupTree const& tree,
                                                   std::vector<proc_t>& procAssignedToEachBlock) const
    {
        // Root node of tree holds total fluid sites
        auto total_sites = tree.levels[0].sites_per_node[0];

        // Those blocks which contain at least one fluid site are here
        auto& nonsolid_block_fluid_site_counts = tree.levels[tree.n_levels].sites_per_node;
        auto const n_nonsolid = nonsolid_block_fluid_site_counts.size();

        if (n_nonsolid < comm_size)
            throw (Exception() << "More MPI processes than blocks - ParMETIS will be unhappy.");
        std::vector<U64> cumulative_fluid_sites(n_nonsolid + 1);
        cumulative_fluid_sites[0] = 0;
        std::inclusive_scan(nonsolid_block_fluid_site_counts.begin(), nonsolid_block_fluid_site_counts.end(),
                            ++cumulative_fluid_sites.begin());

        if (cumulative_fluid_sites.back() != total_sites)
            throw (Exception() << "Octree is inconsistent");

        // Going to divide the blocks amongst the ranks. Start with them all assigned to rank 0
        std::vector<int> rank_for_block(n_nonsolid, 0);
        assign_range(cumulative_fluid_sites.begin(), --cumulative_fluid_sites.end(),
                     rank_for_block.begin(), rank_for_block.end(), comm_size);

        // We need to return data organised in GMY file order
        // Initialise output to SOLID, we will overwrite the non-solid below
        std::fill(procAssignedToEachBlock.begin(), procAssignedToEachBlock.end(), SITE_OR_BLOCK_SOLID);

        // We need to know the octree ID to map to 3D coords.
        auto& block_oct_ids = tree.levels[tree.n_levels].node_ids;
        auto& dims = geometry.GetBlockDimensions();
        auto const strides = BlockLocation{dims[1]*dims[2], dims[2], 1};

        for (std::size_t i = 0; i < n_nonsolid; ++i) {
            auto ijk = octree::oct_to_ijk(block_oct_ids[i]);
            auto gmy = Dot(ijk, strides);
            procAssignedToEachBlock[gmy] = rank_for_block[i];
        }
        return rank_for_block;
    }

    void BasicDecomposition::Validate(std::vector<proc_t>& procAssignedToEachBlock, net::MpiCommunicator const& communicator) const
    {
        log::Logger::Log<log::Debug, log::OnePerCore>("Validating procForEachBlock");

        std::vector<proc_t> procForEachBlockRecv = communicator.AllReduce(procAssignedToEachBlock,
                                                                          MPI_MAX);

        for (site_t block = 0; block < geometry.GetBlockCount(); ++block)
        {
            if (procAssignedToEachBlock[block] != procForEachBlockRecv[block])
            {
                log::Logger::Log<log::Critical, log::OnePerCore>("At least one other proc thought block %li should be on proc %li but we locally had it as %li",
                                                                 block,
                                                                 procAssignedToEachBlock[block],
                                                                 procForEachBlockRecv[block]);
            }
        }
    }

}
