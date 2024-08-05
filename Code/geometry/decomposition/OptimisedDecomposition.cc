// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "geometry/ParmetisHeader.h"
#include "geometry/decomposition/OptimisedDecomposition.h"
#include "geometry/decomposition/DecompositionWeights.h"
#include "geometry/LookupTree.h"

#include "lb/lattices/D3Q27.h"
#include "log/Logger.h"
#include "net/SparseExchange.h"
#include "util/Iterator.h"
#include "util/span.h"

namespace hemelb::geometry::decomposition
{
    using namespace util;

    OptimisedDecomposition::OptimisedDecomposition(
        reporting::Timers& timers,
        net::MpiCommunicator c,
        const GmyReadResult& geometry,
        const lb::LatticeInfo& latticeInfo
    ) : timers(timers), comms(std::move(c)), geometry(geometry),
        tree(geometry.block_store->GetTree()), latticeInfo(latticeInfo),
        procForBlockOct(geometry.block_store->GetBlockOwnerRank()),
        fluidSitesPerBlockOct(tree.levels.back().sites_per_node)
    {
        timers.InitialGeometryRead().Start(); //overall dbg timing

        // Calculate the site distribution and validate if appropriate.
        PopulateSiteDistribution();

        if constexpr (build_info::VALIDATE_GEOMETRY) {
          ValidateVertexDistribution();
          ValidateFirstSiteIndexOnEachBlock();
        }

        // Populate the adjacency data arrays (for ParMetis) and validate if appropriate
        idx_t localVertexCount = vtxDistribn[comms.Rank() + 1] - vtxDistribn[comms.Rank()];

        PopulateAdjacencyData(localVertexCount);

        if constexpr (build_info::VALIDATE_GEOMETRY) {
          ValidateAdjacencyData(localVertexCount);
        }

        log::Logger::Log<log::Trace, log::OnePerCore>("Adj length %i", localAdjacencies.size());

        timers.InitialGeometryRead().Stop();

        // Call parmetis.
        timers.parmetis().Start();
        log::Logger::Log<log::Debug, log::OnePerCore>("Making the call to Parmetis");
        CallParmetis(localVertexCount);
        timers.parmetis().Stop();
        log::Logger::Log<log::Debug, log::OnePerCore>("Parmetis has finished.");

        // Now each process knows which rank all its sites belong
        // on. Tell the destination ranks which sites they need.
        timers.PopulateOptimisationMovesList().Start();
        log::Logger::Log<log::Debug, log::OnePerCore>("Getting moves lists for this core.");
        PopulateMovesList();
        log::Logger::Log<log::Debug, log::OnePerCore>("Done getting moves lists for this core");
        timers.PopulateOptimisationMovesList().Stop();
    }

      void OptimisedDecomposition::CallParmetis(idx_t localVertexCount)
      {
        // From the ParMETIS documentation:
        // --------------------------------
        // Processor Pi holds ni consecutive vertices and mi corresponding edges
        //
        // xadj[ni+1] has the cumulative number of adjacencies per vertex (with a leading 0 on each processor)
        // vwgt[ni] has vertex weight coefficients and can be nullptr
        // adjncy[mi] has the adjacent vertices for each edge (using a global index, starting at 0)
        // adjwgt[mi] has edge weights and can be nullptr
        // vtxdist[P+1] has an identical array of the number of the vertices on each processor, cumulatively.
        //           So Pi has vertices from vtxdist[i] to vtxdist[i+1]-1
        // wgtflag* is 0 with no weights (1 on edges, 2 on vertices, 3 on edges & vertices)
        // numflag* is 0 for C-style numbering (1 for Fortran-style)
        // ncon* is the number of weights on each vertex
        // nparts* is the number of sub-domains (partition domains) desired
        // tpwgts* is the fraction of vertex weight to apply to each sub-domain
        // ubvec* is an array of the imbalance tolerance for each vertex weight
        // options* is an int array of options
        //
        // edgecut[1] will contain the number of edges cut by the partitioning
        // part[ni] will contain the partition vector of the locally-stored vertices
        // comm* is a pointer to the MPI communicator of the processes involved

        // Initialise the partition vector.
        partitionVector = std::vector<idx_t>(localVertexCount, comms.Rank());

        // Weight all vertices evenly.
        //std::vector < idx_t > vertexWeight(localVertexCount, 1);

        // Populate the vertex weight data arrays (for ParMetis) and
        // print out the number of different fluid sites on each core
        PopulateVertexWeightData(localVertexCount);

        // Going to follow ParMETIS docs for naming parameters. Comments to make less cryptic!

        // Desired number of sub domains.
        idx_t ncon = 1; // Number of constraints
        idx_t wgtflag = 2; // Indicate weighting 2 => vertices only
        idx_t numflag = 0; // Numbering scheme 0 => C style from 0
        idx_t nparts = comms.Size();

        // Set the weights of each partition to be even, and to sum to 1.
        std::vector<real_t> tpwgts(nparts, real_t(1.0) / real_t(nparts));
        // A bunch of values ParMetis needs.
        idx_t edgesCut = 0;
        idx_t options[4] = { 0, 0, 0, 0 };
        if constexpr (build_info::VALIDATE_GEOMETRY) {
          // Specify that some options are set and that we should
          // debug everything.
          // Specify that we have set some options
          options[0] = 1;
          // From parmetis.h
          // We get timing info (1)
          // more timing info (2)
          // details of the graph-coarsening process (4)
          // info during graph refinement (8)
          // NOT info on matching (16)
          // info on communication during matching (32)
          // info on remappining (64)
          options[1] = 1 | 2 | 4 | 8 | 32 | 64;
        }
        // Tolerance. 1 => perfect balance, npars => perfect imbalance.
        // Docs recommend 1.05 so we are being quite strict here.
        real_t ubvec = 1.001F;
        MPI_Comm communicator = comms;

        log::Logger::Log<log::Debug, log::OnePerCore>("Calling ParMetis");
        int err = ParMETIS_V3_PartKway(
            vtxDistribn.data(),
            adjacenciesPerVertex.data(),
            localAdjacencies.data(),
            vertexWeights.data(),
            nullptr, // adjwgt, adjacency weights
            &wgtflag,
            &numflag,
            &ncon,
            &nparts,
            tpwgts.data(),
            &ubvec,
            options,
            &edgesCut,
            partitionVector.data(),
            &communicator
        );

        if (err != METIS_OK) {
          throw Exception() << "ParMETIS error";
        }

        log::Logger::Log<log::Debug, log::OnePerCore>("ParMetis returned.");
        if (comms.Rank() == comms.Size() - 1)
        {
          log::Logger::Log<log::Info, log::OnePerCore>("ParMetis cut %d edges.", edgesCut);
          if (edgesCut < 1 && comms.Size() > 2)
          {
            throw Exception()
                << "The decomposition using ParMetis returned an edge cut of 0 even though there are multiple processes. "
                << "This means/implies that ParMETIS cannot properly decompose the system, and no properly load-balanced parallel simulation can be started.";
          }
        }
      }

    void OptimisedDecomposition::PopulateVertexWeightData(idx_t localVertexCount)
    {
        // These counters will be used later on to count the number of each type of vertex site
        std::array<int, 6> siteCounters;
        std::fill(begin(siteCounters), end(siteCounters), 0);

        vertexWeights.resize(localVertexCount);
        idx_t i_wgt = 0;
        // For each block (counting up by lowest site id)...
        for (auto [block_idx, block_rank]: enumerate_with<std::size_t>(procForBlockOct)) {
            if (block_rank != comms.Rank())
                // Only consider sites on this procesr.
                continue;

            auto block_ijk = tree.GetLeafCoords(block_idx);
            auto block_gmy = geometry.GetBlockIdFromBlockCoordinates(block_ijk);

            const BlockReadResult& blockReadResult = geometry.Blocks[block_gmy];

            // util::Vector3D<int> block_coord = BS * block_ijk;
            // ... iterate over sites within the block...
            for (auto [local_site_ijk, m]: IterSitesInBlock(geometry)) {
                // ... only looking at non-solid sites...
                if (blockReadResult.Sites[m].targetProcessor == SITE_OR_BLOCK_SOLID)
                    continue;

                SiteData siteData(blockReadResult.Sites[m]);
                int site_type_i = [&] () {
                    switch (siteData.GetCollisionType()) {
                    case FLUID:
                        return 0;
                    case WALL:
                        return 1;
                    case INLET:
                        return 2;
                    case OUTLET:
                        return 3;
                    case (INLET | WALL):
                        return 4;
                    case (OUTLET | WALL):
                        return 5;
                    default:
                        throw Exception() << "Bad collision type";
                    }
                }();
                ++siteCounters[site_type_i];
                vertexWeights[i_wgt++] = hemelbSiteWeights[site_type_i];
            }
        }

        if (i_wgt != localVertexCount)
          throw Exception() << "Wrong number of vertices: expected " << localVertexCount << " got " << i_wgt;

        int TotalCoreWeight = std::inner_product(begin(siteCounters), end(siteCounters),
                                                 hemelbSiteWeights, 0);
        int TotalSites = std::reduce(begin(siteCounters), end(siteCounters), 0);

        log::Logger::Log<log::Debug, log::OnePerCore>("There are %u Bulk Flow Sites, %u Wall Sites, %u IO Sites, %u WallIO Sites on core %u. Total: %u (Weighted %u Points)",
                                                      siteCounters[0],
                                                      siteCounters[1],
                                                      siteCounters[2] + siteCounters[3],
                                                      siteCounters[4] + siteCounters[5],
                                                      comms.Rank(),
                                                      TotalSites,
                                                      TotalCoreWeight);
    }

    void OptimisedDecomposition::PopulateSiteDistribution()
    {
        // Parmetis needs to know (on all processes) how
        // vertices/sites are distributed across the processes.
        //
        // The `vtxdist` array has P+1 elements.
        //
        // `vtxdist[i]` gives the total number of vertices that are on
        // processes with rank less than `i`.
        //
        // `vtxdist[0]` is always zero
        //
        // `vtxdist[P]` has the total number of vertices
        //
        // `vtxdist[i+1] - vtxdist[i]` gives the number of vertices on
        // process i
        //
        // We also need to know what vertex indices to assign to each
        // site. Do this contiguously for the sites within each block
        // per-process.
        //
        // We will also need to be able to translate between a site's
        // GMY index in a block and it's fluid-only index.

        auto& procForBlockOct = geometry.block_store->GetBlockOwnerRank();
        auto& tree = geometry.block_store->GetTree();
        auto& fluidSitesPerBlockOct = tree.levels[tree.n_levels].sites_per_node;
        auto const NBLOCKS = geometry.block_store->GetBlockCount();

        // First, count the sites per process and assign contiguous
        // IDs to the sites in each block.
        U64 total_sites = 0;
        vtxCountPerProc = std::vector<idx_t>(comms.Size(), 0);
        firstSiteIndexPerBlockOct.resize(NBLOCKS + 1);

        for (site_t block = 0; block < NBLOCKS; ++block) {
            auto block_rank = procForBlockOct[block];
            firstSiteIndexPerBlockOct[block] = total_sites;
            vtxCountPerProc[block_rank] += fluidSitesPerBlockOct[block];
            total_sites += fluidSitesPerBlockOct[block];
        }
        firstSiteIndexPerBlockOct[NBLOCKS] = total_sites;

        // Now do the scan to setup vtxdist
        vtxDistribn = std::vector<idx_t>(comms.Size() + 1);
        vtxDistribn[0] = 0;
        std::inclusive_scan(
            vtxCountPerProc.begin(), vtxCountPerProc.end(),
            vtxDistribn.begin() + 1,
            std::plus<idx_t>()
        );

        // Now, for every block we've got data for, create the mapping
        // from block OCT id to a vector of the GMY local site ids of
        // the fluid sites.
        for (auto [block_idx, block_rank]: enumerate_with<U64>(procForBlockOct)) {
            auto block_ijk = tree.GetLeafCoords(block_idx);
            auto block_gmy = geometry.GetBlockIdFromBlockCoordinates(block_ijk);
            auto& blockReadResult = geometry.Blocks[block_gmy];
            if (blockReadResult.Sites.empty())
                continue;
            std::vector<U16>& block_fluid_site_gmy_ids = gmySiteIdForBlockOct[block_idx];
            for (auto const& [i, s]: util::enumerate_with<U16>(blockReadResult.Sites)) {
                // ... only looking at non-solid sites...
                if (s.targetProcessor != SITE_OR_BLOCK_SOLID) {
                    block_fluid_site_gmy_ids.push_back(i);
                }
            }
        }
    }

    void OptimisedDecomposition::PopulateAdjacencyData(idx_t localVertexCount)
    {
        adjacenciesPerVertex.push_back(0);

        auto const& block_dims = geometry.GetBlockDimensions();
        auto const BS = geometry.GetBlockSize();
        // For range check below
        auto const lo_site = Vector3D<site_t>::Zero();
        auto const hi_site = block_dims.as<site_t>() * BS - Vector3D<site_t>::Ones();

        // For each block (counting up by lowest site id)...
        for (auto [block_idx, block_rank]: enumerate_with<std::size_t>(procForBlockOct)) {
            if (block_rank != comms.Rank())
                // ... considering only the ones which live on this proc...
                continue;

            auto block_ijk = tree.GetLeafCoords(block_idx);
            auto block_gmy = geometry.GetBlockIdFromBlockCoordinates(block_ijk);

            auto& blockReadResult = geometry.Blocks[block_gmy];
            for (auto [local_site_ijk, m]: IterSitesInBlock(geometry)) {
                // ... only looking at non-solid sites...
                if (blockReadResult.Sites[m].targetProcessor == SITE_OR_BLOCK_SOLID)
                    continue;

                auto global_site_ijk = block_ijk.as<site_t>() * BS + local_site_ijk;
                // ... for each lattice direction (recalling vector[0] == {0,0,0}) ...
                for (unsigned int l = 1; l < latticeInfo.GetNumVectors(); l++) {
                    // ... which leads to a valid neighbouring site...
                    auto neigh_ijk = global_site_ijk + latticeInfo.GetVector(l);
                    if (!neigh_ijk.IsInRange(lo_site, hi_site))
                        continue;

                    // ... (that is actually being simulated and not a solid)...
                    auto neigh_block = Vec16(neigh_ijk / BS);
                    auto neigh_idx = tree.GetPath(neigh_block).leaf();
                    // Recall that the octree will return a path
                    // with "no child" for all levels where the
                    // requested node doesn't exist.
                    if (neigh_idx == octree::Level::NC)
                        continue;

                    auto [ni, nj, nk] = neigh_ijk % BS;
                    U16 neighbourSiteId = geometry.GetSiteIdFromSiteCoordinates(ni, nj, nk);

                    // This is the GMY local site IDs for only the
                    // fluid sites on the neighbour block
                    auto const& neigh_gmy_sites = gmySiteIdForBlockOct[neigh_idx];
                    auto lb = std::lower_bound(neigh_gmy_sites.begin(), neigh_gmy_sites.end(), neighbourSiteId);
                    // lb is either end or the first elem greater than or equal to what we want
                    // end => SOLID; greater than or equal => SOLID
                    if (lb == neigh_gmy_sites.end() || *lb != neighbourSiteId)
                        continue;
                    // have equal
                    auto local_contig_idx = std::distance(neigh_gmy_sites.begin(), lb);
                    U64 neighGlobalSiteId = firstSiteIndexPerBlockOct[neigh_idx] + local_contig_idx;

                    // then add this to the list of adjacencies.
                    localAdjacencies.push_back(idx_t(neighGlobalSiteId));
                }

                // The cumulative count of adjacencies for this vertex is equal to the total
                // number of adjacencies we've entered.
                // NOTE: The prefix operator is correct here because
                // the array has a leading 0 not relating to any site.
                adjacenciesPerVertex.push_back(localAdjacencies.size());
            }
        }
    }

    // We are returning a map (keyed by the process rank) of where
    // each site that starts on this process should end up. Every site
    // is here and the values of the map (vectors) will be sorted,
    // first by block ID and then by site ID (note this for use in the
    // geometry reader!)
    auto OptimisedDecomposition::CompileMoveData() -> MovesMap
    {
        // Key is the destination rank, value is the list of [block, site id] pairs
        std::map<int, std::vector<std::array<U64,2>>> movesToRank;

        std::vector<idx_t> moveData;
        auto const myLowest = vtxDistribn[comms.Rank()];
        auto const mySiteCount = vtxDistribn[comms.Rank() + 1] - myLowest;
        // For each local fluid site...
        // (NB: since the vertices/sites are sorted by block and then
        // site, the results are also so sorted.)
        for (idx_t ii = 0; ii < mySiteCount; ++ii)
        {
            // Where should it be? Catalogue even sites staying here for simplicity
            auto dest_rank = partitionVector[ii];

            idx_t global_site_id = myLowest + ii;
            // Find out which block it's on.

            // Upper bound gives the iterator to the first one after
            // (recalling that we have an extra entry with the total
            // number of the fluid sites at the end).
            auto after = std::upper_bound(
                firstSiteIndexPerBlockOct.begin(), firstSiteIndexPerBlockOct.end(),
                global_site_id
            );
            U64 const block_idx = std::distance(firstSiteIndexPerBlockOct.begin(), after) - 1;

            if constexpr (build_info::VALIDATE_GEOMETRY) {
                // Check the block id is correct
              if (
                  global_site_id  < firstSiteIndexPerBlockOct[block_idx] ||
                  global_site_id  >= firstSiteIndexPerBlockOct[block_idx+1]) {
                  log::Logger::Log<log::Critical, log::OnePerCore>(
                      "Found block %i for site %i but sites on this block start at number %i and finish before %i",
                      block_idx,
                      global_site_id ,
                      firstSiteIndexPerBlockOct[block_idx],
                      firstSiteIndexPerBlockOct[block_idx + 1]
                  );
              }
            }

            // ... and find its site id within that block. Start by working out how many fluid sites
            // we have to pass before we arrive at the fluid site we're after...
            site_t block_site_id = global_site_id - firstSiteIndexPerBlockOct[block_idx];
            auto siteIndexGmy = gmySiteIdForBlockOct[block_idx][block_site_id];

            // Add the block, site and destination rank to our move list.
            movesToRank[dest_rank].push_back({block_idx, siteIndexGmy});
        }

        return movesToRank;
      }

    // Organise out the moves: those sites stay, leaving and arriving.
    void OptimisedDecomposition::PopulateMovesList() {
        auto rank = comms.Rank();

        // Work out where *all* sites should be
        leaving = CompileMoveData();

        // Separate out those staying
        staying = std::move(leaving[rank]);
        leaving.erase(rank);

        // Spread the move data around
        log::Logger::Log<log::Debug, log::OnePerCore>("Starting to spread move data");

        net::sparse_exchange<SiteDesc> xchg(comms, 444);
        for (auto const& [dest, moves]: leaving) {
            xchg.send(to_span(moves), dest);
        }

        xchg.receive(
            [&](int src, int count) {
                auto& rbuf = arriving[src];
                rbuf.resize(count);
                return rbuf.data();
            },
            [&](int src, auto* buf) {
                // no-op
            }
        );
    }

    void OptimisedDecomposition::ValidateVertexDistribution() {
        log::Logger::Log<log::Debug, log::OnePerCore>("Validating the vertex distribution.");
        // vtxDistribn should be the same on all cores.
        std::vector<idx_t> vtxDistribnRecv = comms.AllReduce(vtxDistribn, MPI_MIN);

        for (proc_t rank = 0; rank < comms.Size() + 1; ++rank) {
            if (vtxDistribn[rank] != vtxDistribnRecv[rank]) {
                log::Logger::Log<log::Critical, log::OnePerCore>(
                    "vtxDistribn[%i] was %li but at least one other core had it as %li.",
                    rank,
                    vtxDistribn[rank],
                    vtxDistribnRecv[rank]
                );
            }
        }
    }

    void OptimisedDecomposition::ValidateAdjacencyData(idx_t localVertexCount) {
        // If we're using debugging logs, check that the arguments are consistent across all cores.
        // To verify: vtxDistribn, adjacenciesPerVertex, adjacencies
        log::Logger::Log<log::Debug, log::OnePerCore>("Validating the graph adjacency structure");
        // Create an array of lists to store all of this node's adjacencies, arranged by the
        // proc the adjacent vertex is on.
        using edge = std::array<idx_t, 2>;
        std::map<int, std::vector<edge>> adjByNeighProc;
        auto check_edge = [this] (idx_t v, idx_t w) {
          idx_t w_local = w - vtxDistribn[comms.Rank()];
          bool found = false;
          for (
              idx_t edge_i = adjacenciesPerVertex[w_local];
              edge_i < adjacenciesPerVertex[w_local + 1];
              ++edge_i
          ) {
            if (localAdjacencies[edge_i] == v) {
              found = true;
              break;
            }
          }
          if (!found)
            log::Logger::Log<log::Critical, log::OnePerCore>(
                "Could not find reverse for edge from vertex %li to %li",
                v, w
            );
        };

        // The adjacency data should correspond across all cores.
        for (idx_t index = 0; index < localVertexCount; ++index) {
            idx_t vertex = vtxDistribn[comms.Rank()] + index;
            // Iterate over each adjacency (of each vertex).
            for (
                idx_t adj_idx = adjacenciesPerVertex[index];
                adj_idx < adjacenciesPerVertex[index + 1];
                ++adj_idx
            ) {
                idx_t adjacentVertex = localAdjacencies[adj_idx];
                // Find the proc of the neighbouring vertex.
                // This iterator should point to the next process's start pos
                auto ub = std::upper_bound(vtxDistribn.begin(), vtxDistribn.end(), adjacentVertex);
                if (ub == vtxDistribn.end()) {
                    log::Logger::Log<log::Critical, log::OnePerCore>("The vertex %li has a neighbour %li which doesn\'t appear to live on any processor.",
                                                                     vertex,
                                                                     adjacentVertex);
                    continue;
                }
                int adj_proc = std::distance(vtxDistribn.begin(), ub) - 1;

                // Most edges should be local - check these in place to speed it up.
                if (adj_proc == comms.Rank()) {
                  check_edge(vertex, adjacentVertex);
                } else {
                  // Store the edge for later checking
                  adjByNeighProc[adj_proc].push_back({adjacentVertex, vertex});
                }
            }
        }

        // Now sort the cross-process edges
        for (auto& [_, edges]: adjByNeighProc) {
            std::sort(
                edges.begin(), edges.end(), [](edge const& a, edge const& b) {
                  return a[0] < b[0];
                }
            );
        }

        std::vector<std::size_t> neigh_edge_counts(comms.Size());
        for (int root = 0; root < comms.Size(); ++root) {
          // Send any edges we've got to the current root process
          auto maybe_edges = adjByNeighProc.find(root);
          std::size_t const nedges = (maybe_edges == adjByNeighProc.end()) ?
            0 // no edges to that process
            : maybe_edges->second.size();

          net::MpiCall{MPI_Gather}(&nedges, 1, net::MpiDataType<std::size_t>(), neigh_edge_counts.data(), 1, net::MpiDataType<std::size_t>(), root, comms);

          if (root == comms.Rank()) {
            // Root receives and checks non-zeros
            std::vector<edge> neigh_edges;
            for (int src = 0; src < comms.Size(); ++src) {
              if (neigh_edge_counts[src]) {
                neigh_edges.resize(neigh_edge_counts[src]);
                comms.Receive(neigh_edges, src, 43);
                // Check!
                for (auto [w, v]: neigh_edges) {
                  check_edge(v, w);
                }
              }
            }
          } else {
            if (nedges)
              comms.Send(maybe_edges->second, root, 43);
          }
        }
    }

    void OptimisedDecomposition::ValidateFirstSiteIndexOnEachBlock() {
        // Check that values are
        // a) the same across all processes
        // b) increase from zero to the the total number of fluid
        // sites by at least one per block and at most by sites per
        // block

        log::Logger::Log<log::Debug, log::OnePerCore>("Validating the firstSiteIndexPerBlock values.");
        auto const n_non_solid = tree.levels[tree.n_levels].sites_per_node.size();
        auto const total_fluid = tree.levels[0].sites_per_node[0];

        HASSERT(firstSiteIndexPerBlockOct.size() == n_non_solid + 1);
        HASSERT(firstSiteIndexPerBlockOct[0] == 0);
        for (auto i = 0UL; i < n_non_solid; ++i) {
            auto d = firstSiteIndexPerBlockOct[i+1] - firstSiteIndexPerBlockOct[i];
            HASSERT(d > 0);
            HASSERT(d <= geometry.GetSitesPerBlock());
        }
        HASSERT(firstSiteIndexPerBlockOct[n_non_solid] == total_fluid);

        auto firstSiteIndexPerBlockRecv = comms.AllReduce(firstSiteIndexPerBlockOct, MPI_MAX);

        for (auto i = 0UL; i < n_non_solid; ++i) {
            if (firstSiteIndexPerBlockOct[i] != firstSiteIndexPerBlockRecv[i]) {
                log::Logger::Log<log::Critical, log::OnePerCore>(
                    "This process had the first site index on block %li as %li but at least one other core had it as %li.",
                    i,
                    firstSiteIndexPerBlockOct[i],
                    firstSiteIndexPerBlockRecv[i]
                );
            }
        }
    }
}
