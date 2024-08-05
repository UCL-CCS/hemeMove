// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "log/Logger.h"
#include "geometry/BlockTraverser.h"
#include "geometry/Domain.h"
#include "geometry/GmyReadResult.h"
#include "geometry/neighbouring/NeighbouringDomain.h"
#include "geometry/LookupTree.h"
#include "net/IOCommunicator.h"
#include "net/net.h"
#include "reporting/Dict.h"

namespace hemelb::geometry
{
        Domain::Domain(const lb::LatticeInfo& latticeInfo,
                       const net::IOCommunicator& comms_) :
                latticeInfo(latticeInfo),
                shared_site_type_indices(comms_, 0),
                neighbouringData(new neighbouring::NeighbouringDomain(latticeInfo)), comms(comms_)
        {
        }

        Domain::Domain(const lb::LatticeInfo& latticeInfo,
                       GmyReadResult& readResult, const net::IOCommunicator& comms_) :
                latticeInfo(latticeInfo), shared_site_type_indices(comms_, 0),
                neighbouringData(new neighbouring::NeighbouringDomain(latticeInfo)),
                rank_for_site_store(std::move(readResult.block_store)),
                comms(comms_)
        {
            SetBasicDetails(readResult.GetBlockDimensions(), readResult.GetBlockSize());

            ProcessReadSites(readResult);
            // if debugging then output beliefs regarding geometry and neighbour list
            if (log::Logger::ShouldDisplay<log::Trace>())
            {
                proc_t localRank = comms.Rank();
                for (auto& neighbouringProc: neighbouringProcs)
                {
                    log::Logger::Log<log::Trace, log::OnePerCore>("domain_type: Rank %i thinks that rank %i is a neighbour with %i shared edges\n",
                                                                  localRank,
                                                                  neighbouringProc.Rank,
                                                                  neighbouringProc.SharedDistributionCount);
                }
            }
            CollectFluidSiteDistribution();
            CollectGlobalSiteExtrema();
            InitialiseNeighbourLookups();
        }

    std::size_t Domain::GetBlockOctIndexFromBlockCoords(const util::Vector3D<std::uint16_t> &blockCoords) const {
        return rank_for_site_store->GetTree().GetLeaf(blockCoords).index();
    }

        // Need d'tor in its own TU to avoid include of LookupTree
        Domain::~Domain() noexcept = default;

        void Domain::SetBasicDetails(Vec16 blocksIn, U16 blockSizeIn)
        {
            blockCounts = blocksIn;
            blockSize = blockSizeIn;
            sites = blocksIn.as<site_t>() * blockSize;
            sitesPerBlockVolumeUnit = blockSize * blockSize * blockSize;
            blockCount = blockCounts.x() * blockCounts.y() * blockCounts.z();
        }

    void Domain::ProcessReadSites(const GmyReadResult & readResult)
    {
        log::Logger::Log<log::Info, log::Singleton>("Processing sites assigned to each MPI process");
	const auto NO_NORMAL = util::Vector3D<float>(NO_VALUE);
        const auto max_site_index = sites - util::Vector3D<site_t>::Ones();
        blocks.resize(rank_for_site_store->GetBlockCount());

        totalSharedFs = 0;

        ReadDataBundle domainEdge;
        ReadDataBundle midDomain;

        proc_t localRank = comms.Rank();

        // We need to know the coordinates of all the domain edge
        // sites below so we can loop over them to figure out which
        // processes they live on once the RMA data structures are
        // initialised.
        std::vector<util::Vector3D<site_t>> edge_sites;

        // This closure returns which rank a site (given by global
        // 3D coordinate) lives on, according to the
        // GmyReadResult. Since that is only guaranteed to know
        // the rank of sites which are assigned to this rank, it
        // may well return UNKNOWN_PROCESS
        auto read_result_rank_for_site = [&max_site_index, &readResult](
            util::Vector3D<site_t> const& global_site_coords,
            Vec16& block_coords, site_t& site_gmy_idx
        ) {
            // If outside the bounding box it is definitely non-fluid
            if (!global_site_coords.IsInRange(
                util::Vector3D<site_t>::Zero(),
                max_site_index
            ))
                return SITE_OR_BLOCK_SOLID;

          // ... (that is actually being simulated and not a solid)...
          block_coords = (global_site_coords / readResult.GetBlockSize()).as<U16>();
          site_t block_gmy_idx = readResult.GetBlockIdFromBlockCoordinates(block_coords.x(),
                                                                           block_coords.y(),
                                                                           block_coords.z());

          // Move on if the neighbour is in a block of solids
          // in which case the block will contain zero sites
          // Or on if the neighbour site is solid
          // in which case the targetProcessor is SITE_OR_BLOCK_SOLID
          // Or the neighbour is also on this processor
          // in which case the targetProcessor is localRank
          if (readResult.Blocks[block_gmy_idx].Sites.empty())
              return SITE_OR_BLOCK_SOLID;

          auto site_local_coords = global_site_coords % readResult.GetBlockSize();
          site_gmy_idx = readResult.GetSiteIdFromSiteCoordinates(site_local_coords.x(),
                                                                 site_local_coords.y(),
                                                                 site_local_coords.z());
          return readResult.Blocks[block_gmy_idx].Sites[site_gmy_idx].targetProcessor;
        };

        for (auto leaf: rank_for_site_store->GetTree().IterLeaves()) {
            auto block_ijk = leaf.coords();
            site_t blockGmyIdx = GetBlockGmyIdxFromBlockCoords(block_ijk);
            auto const& blockReadIn = readResult.Blocks[blockGmyIdx];

            if (blockReadIn.Sites.empty())
                continue;

            auto blockOctIdx = leaf.index();
            auto lowest_site_in_block = block_ijk.as<site_t>() * GetBlockSize();

            // Iterate over all sites within the current block.
            for (
                auto siteTraverser = SiteTraverser(*this);
                siteTraverser.CurrentLocationValid();
                siteTraverser.TraverseOne()
             ) {
                site_t localSiteId = siteTraverser.GetCurrentIndex();
                // Create block if required
                if (blocks[blockOctIdx].IsEmpty())
                    blocks[blockOctIdx] = Block(GetSitesPerBlockVolumeUnit());

		auto const& siteReadIn = blockReadIn.Sites[localSiteId];
                auto assignedRank = siteReadIn.targetProcessor;
                blocks[blockOctIdx].SetProcessorRankForSite(localSiteId, assignedRank);

                // If the site is not on this processor, continue.
                if (assignedRank != localRank)
                    continue;

                bool isMidDomainSite = true;
                util::Vector3D<site_t> siteGlobalCoords = lowest_site_in_block + siteTraverser.GetCurrentLocation();
                // Iterate over all non-zero direction vectors.
                for (unsigned int l = 1; l < latticeInfo.GetNumVectors(); l++) {
                    // Find the neighbour site coords in this direction.
                    util::Vector3D<site_t> neighbourGlobalCoords = siteGlobalCoords
                      + latticeInfo.GetVector(l).as<site_t>();
                    Vec16 neighbourBlock;
                    site_t neighbourSiteId;
                    site_t neighbourProc = read_result_rank_for_site(
                        neighbourGlobalCoords, neighbourBlock, neighbourSiteId
                    );
                    if (neighbourProc == SITE_OR_BLOCK_SOLID || localRank == neighbourProc)
                        continue;
                    isMidDomainSite = false;
                    totalSharedFs++;
                }

                if (!isMidDomainSite)
                    edge_sites.push_back(siteGlobalCoords);

                // Set the collision type data. map_block site data is renumbered according to
                // fluid site numbers within a particular collision type.
                SiteData siteData(siteReadIn);
                int l = -1;
                switch (siteData.GetCollisionType()) {
                case FLUID:
                    l = 0;
                    break;
                case WALL:
                    l = 1;
                    break;
                case INLET:
                    l = 2;
                    break;
                case OUTLET:
                    l = 3;
                    break;
                case (INLET | WALL):
                    l = 4;
                    break;
                case (OUTLET | WALL):
                    l = 5;
                    break;
                }

                auto& normal = siteReadIn.wallNormalAvailable ? siteReadIn.wallNormal : NO_NORMAL;

                auto push_site_onto_bundle = [&](ReadDataBundle& b) {
                    b.blockNumbers[l].push_back(blockOctIdx);
                    b.siteNumbers[l].push_back(localSiteId);
                    b.siteData[l].push_back(siteData);
                    b.wallNormals[l].push_back(normal);
                    for (Direction direction = 1; direction < latticeInfo.GetNumVectors(); direction++) {
                        b.wallDistance[l].push_back(siteReadIn.links[direction - 1].distanceToIntersection);
                    }
                };
                if (isMidDomainSite) {
                    push_site_onto_bundle(midDomain);
                } else {
                    push_site_onto_bundle(domainEdge);
                }
            }
        }

        PopulateWithReadData(midDomain, domainEdge);

        // We now have the distributed store setup, so can find which
        // process any site lives one. Set up the neighbouring
        // processes for the edge sites.
        for (auto const& site_global_coords: edge_sites) {
            for (unsigned int l = 1; l < latticeInfo.GetNumVectors(); l++) {
                // Find the neighbour site co-ords in this direction.
                util::Vector3D<site_t> neigh_global_coords = site_global_coords + latticeInfo.GetVector(l).as<site_t>();
                Vec16 neigh_block;
                site_t neigh_site_id;
                site_t rrProc = read_result_rank_for_site(
                    neigh_global_coords, neigh_block, neigh_site_id
                );

                if (rrProc == SITE_OR_BLOCK_SOLID || rrProc == localRank)
                    continue;

                auto [neighbourProc, remoteSiteIdx] = rank_for_site_store->GetSiteData(
                    neigh_block, neigh_site_id
                );
                auto neighProcWithSite = std::find_if(
                    neighbouringProcs.begin(), neighbouringProcs.end(),
                    [&](NeighbouringProcessor const& np) {
                      return np.Rank == neighbourProc;
                    }
                );

                if (neighProcWithSite == neighbouringProcs.end()) {
                    // We didn't find a neighbour-proc with the
                    // neighbour-site on it, so we need a new
                    // neighbouring processor.

                    // Store rank of neighbour in >neigh_proc[neigh_procs]
                    NeighbouringProcessor lNewNeighbour;
                    lNewNeighbour.SharedDistributionCount = 1;
                    lNewNeighbour.Rank = neighbourProc;
                    neighbouringProcs.push_back(lNewNeighbour);
                } else {
                    // Did find it, increment the shared count
                    ++neighProcWithSite->SharedDistributionCount;
                }
            }
        }
    }

    void Domain::PopulateWithReadData(const ReadDataBundle &midDomain, const ReadDataBundle &domainEdge)
    {
        auto write_my_sites = rank_for_site_store->begin_writes();

        log::Logger::Log<log::Info, log::Singleton>("Assigning local indices to sites and associated data");

        // Data about local sites.
        SiteRankIndex rank_index = {comms.Rank(), 0};
        // Track index of changes in type and store them in the counts array
        auto ranges = shared_site_type_indices.Span();
        auto current_count_iter = begin(ranges);
        // Merge local sites into one big array. This lambda does either all the
        // midDomain or all the domainEdge sites, putting them in order of their
        // collision type (bulk/wall/inlet/etc).
        auto op = [this, &rank_index, &write_my_sites, &current_count_iter](ReadDataBundle const& b) {
            auto const NV = latticeInfo.GetNumVectors();
            auto& localFluidSites = rank_index[1];
            for (unsigned collisionType = 0; collisionType < COLLISION_TYPES; collisionType++) {
                *current_count_iter = localFluidSites;

                for (unsigned i = 0; i < b.blockNumbers[collisionType].size(); ++i) {
                    siteData.push_back(b.siteData[collisionType][i]);
                    wallNormalAtSite.emplace_back(b.wallNormals[collisionType][i]);
                    for (Direction dir = 1; dir < NV; dir++) {
                        distanceToWall.push_back(
                                b.wallDistance[collisionType][i * (NV - 1) + dir - 1]
                        );
                    }
                    site_t blockId = b.blockNumbers[collisionType][i];
                    site_t siteId = b.siteNumbers[collisionType][i];
                    blocks[blockId].SetLocalContiguousIndexForSite(siteId, localFluidSites);
                    globalSiteCoords.push_back(GetGlobalCoords(blockId, GetSiteCoordsFromSiteId(siteId)));
                    write_my_sites(blockId)(siteId) = rank_index;
                    localFluidSites++;
                }

                // Note for all iterations except the last we will overwrite this, but nice and clear.
                ++current_count_iter;
                *current_count_iter = localFluidSites;
            }
        };
        // Now apply it. Order of these is important!
        op(midDomain);
        op(domainEdge);
    }

        void Domain::CollectFluidSiteDistribution()
        {
            log::Logger::Log<log::Debug, log::Singleton>("Gathering site counts.");
            fluidSitesOnEachProcessor = comms.AllGather(GetLocalFluidSiteCount());
            totalFluidSites = std::reduce(
                    fluidSitesOnEachProcessor.begin(), fluidSitesOnEachProcessor.end(),
                    site_t(0), std::plus<>{});
        }

        void Domain::CollectGlobalSiteExtrema()
        {
            log::Logger::Log<log::Debug, log::Singleton>("Gathering bounds.");
            auto localMins = util::Vector3D<site_t>::Largest();
            auto localMaxes = util::Vector3D<site_t>::Zero();

            for (auto leaf: rank_for_site_store->GetTree().IterLeaves()) {
                auto blockOctIdx = leaf.index();

                auto const& block = blocks[blockOctIdx];
                if (block.IsEmpty())
                    // we may not have read a block on this rank
                    continue;
                auto block_ijk = leaf.coords();
                auto lowest_site_in_block = block_ijk.as<site_t>() * GetBlockSize();
                for (auto siteSet = SiteTraverser(*this);
                     siteSet.CurrentLocationValid(); siteSet.TraverseOne())
                {
                    if (block.GetProcessorRankForSite(siteSet.GetCurrentIndex()) == comms.Rank())
                    {
                        auto globalCoords = lowest_site_in_block + siteSet.GetCurrentLocation();

                        localMins.UpdatePointwiseMin(globalCoords);
                        localMaxes.UpdatePointwiseMax(globalCoords);
                    }
                }

            }

            comms.AllReduceInPlace(std::span(localMins.begin(), localMins.end()), MPI_MIN);
            comms.AllReduceInPlace(std::span(localMaxes.begin(), localMaxes.end()), MPI_MAX);

            globalSiteMins = localMins;
            globalSiteMaxes = localMaxes;
        }

        void Domain::InitialiseNeighbourLookups()
        {
            log::Logger::Log<log::Info, log::Singleton>("Initialising neighbour lookups");
            // Allocate the index in which to put the distribution functions received from the other
            // process.
            //auto sharedDistributionLocationForEachProc = std::vector<std::vector<site_t> >(comms.Size());
            site_t totalSharedDistributionsSoFar = 0;
            // Set the remaining neighbouring processor data.
            for (auto & neighbouringProc : neighbouringProcs)
            {
                // Pointing to a few things, but not setting any variables.
                // FirstSharedF points to start of shared_fs.
                neighbouringProc.FirstSharedDistribution = GetLocalFluidSiteCount()
                                                                         * latticeInfo.GetNumVectors() + 1 + totalSharedDistributionsSoFar;
                totalSharedDistributionsSoFar += neighbouringProc.SharedDistributionCount;
            }
            auto sharedDistributionLocationForEachProc = InitialiseNeighbourLookup();
            InitialisePointToPointComms(sharedDistributionLocationForEachProc);
            InitialiseReceiveLookup(sharedDistributionLocationForEachProc);
        }

        auto Domain::InitialiseNeighbourLookup() -> proc2neighdata
        {
            proc2neighdata ans;
            const proc_t localRank = comms.Rank();
            neighbourIndices.resize(latticeInfo.GetNumVectors() * GetLocalFluidSiteCount());
            for (auto leaf: rank_for_site_store->GetTree().IterLeaves()) {
                auto const& map_block_p = blocks[leaf.index()];
                if (map_block_p.IsEmpty())
                    continue;

                auto lowest_site_in_block = leaf.coords().as<site_t>() * GetBlockSize();
                for (auto siteTraverser = SiteTraverser(*this);
                     siteTraverser.CurrentLocationValid(); siteTraverser.TraverseOne())
                {
                    if (localRank != map_block_p.GetProcessorRankForSite(siteTraverser.GetCurrentIndex()))
                        continue;

                    // Get site data, which is the number of the fluid site on this proc..
                    site_t localIndex =
                            map_block_p.GetLocalContiguousIndexForSite(siteTraverser.GetCurrentIndex());

                    auto currentLocationCoords = lowest_site_in_block + siteTraverser.GetCurrentLocation();
                    // Set neighbour location for the distribution component at the centre of
                    // this site.
                    SetNeighbourLocation(localIndex, 0, localIndex * latticeInfo.GetNumVectors() + 0);
                    for (Direction direction = 1; direction < latticeInfo.GetNumVectors(); direction++)
                    {
                        // Work out positions of neighbours.
                        auto neighbourCoords = currentLocationCoords + latticeInfo.GetVector(direction).as<site_t>();
                        if (!IsValidLatticeSite(neighbourCoords))
                        {
                            // Set the neighbour location to the rubbish site.
                            SetNeighbourLocation(localIndex,
                                                 direction,
                                                 GetLocalFluidSiteCount() * latticeInfo.GetNumVectors());
                            continue;
                        }
                        // Get the id of the processor which the neighbouring site lies on.
                        const proc_t proc_id_p = GetProcIdFromGlobalCoords(neighbourCoords);
                        if (proc_id_p == SITE_OR_BLOCK_SOLID)
                        {
                            // initialize f_id to the rubbish site.
                            SetNeighbourLocation(localIndex,
                                                 direction,
                                                 GetLocalFluidSiteCount() * latticeInfo.GetNumVectors());
                            continue;
                        }
                        else
                            // If on the same proc, set f_id of the
                            // current site and direction to the
                            // site and direction that it sends to.
                            // If we check convergence, the data for
                            // each site is split into that for the
                            // current and previous cycles.
                        if (localRank == proc_id_p)
                        {
                            // Pointer to the neighbour.
                            site_t contigSiteId = GetContiguousSiteId(neighbourCoords);
                            SetNeighbourLocation(localIndex,
                                                 direction,
                                                 contigSiteId * latticeInfo.GetNumVectors() + direction);
                            continue;
                        }
                        else
                        {
                            // This stores some coordinates.  We
                            // still need to know the site number.
                            // neigh_proc[ n ].f_data is now
                            // set as well, since this points to
                            // f_data.  Every process has data for
                            // its neighbours which say which sites
                            // on this process are shared with the
                            // neighbour.
                            ans[proc_id_p].emplace_back(currentLocationCoords, direction);
                        }
                    }
                }

            }
            return ans;
        }

        void Domain::InitialisePointToPointComms(proc2neighdata& sharedFLocationForEachProc)
        {
            proc_t localRank = comms.Rank();
            // point-to-point communications are performed to match data to be
            // sent to/receive from different partitions; in this way, the
            // communication of the locations of the interface-dependent fluid
            // sites and the identifiers of the distribution functions which
            // propagate to different partitions is avoided (only their values
            // will be communicated). It's here!
            // Allocate the request variable.
            std::vector<net::MpiRequest> reqs(neighbouringProcs.size());
            int i_req = 0;
            for (auto& neighbouringProc : neighbouringProcs)
            {
                // One way send receive.  The lower numbered netTop->ProcessorCount send and the higher numbered ones receive.
                // It seems that, for each pair of processors, the lower numbered one ends up with its own
                // edge sites and directions stored and the higher numbered one ends up with those on the
                // other processor.
                if (neighbouringProc.Rank > localRank)
                {
                    // We know that the elements are contiguous from asserts
                    // in Domain.h about size and alignment of point_direction.
                    // Using a template as want the const/mutable variants.
                    auto const& to_send = sharedFLocationForEachProc.at(neighbouringProc.Rank);

                    reqs[i_req] = comms.Issend(
                            std::span<site_t const>(&to_send[0].first[0], to_send.size() * 4),
                            neighbouringProc.Rank
                            );
                } else {
                    auto& dest = sharedFLocationForEachProc[neighbouringProc.Rank];
                    dest.resize(neighbouringProc.SharedDistributionCount);
                    reqs[i_req] = comms.Irecv(
                            std::span<site_t>(&dest[0].first[0], 4*dest.size()),
                            neighbouringProc.Rank
                    );
                }
                i_req += 1;
            }
            net::MpiRequest::Waitall(reqs);
        }

        void Domain::InitialiseReceiveLookup(proc2neighdata const& sharedFLocationForEachProc)
        {
            proc_t localRank = comms.Rank();
            streamingIndicesForReceivedDistributions.resize(totalSharedFs);
            site_t f_count = GetLocalFluidSiteCount() * latticeInfo.GetNumVectors();
            site_t sharedSitesSeen = 0;
            for (auto& neighbouringProc: neighbouringProcs) {
                for (site_t sharedDistributionId = 0;
                     sharedDistributionId < neighbouringProc.SharedDistributionCount; sharedDistributionId++)
                {
                    // Get coordinates and direction of the distribution function to be sent to another process.
                    auto [location, l] = sharedFLocationForEachProc.at(neighbouringProc.Rank)[sharedDistributionId];
                    // Correct so that each process has the correct coordinates.
                    if (neighbouringProc.Rank < localRank)
                    {
                        location += latticeInfo.GetVector(l);
                        l = latticeInfo.GetInverseIndex(l);
                    }
                    // Get the fluid site number of site that will send data to another process.
                    site_t contigSiteId = GetContiguousSiteId(location);
                    // Set f_id to the element in the send buffer that we put the updated
                    // distribution functions in.
                    SetNeighbourLocation(contigSiteId, (unsigned int) ( (l)), ++f_count);
                    // Set the place where we put the received distribution functions, which is
                    // f_new[number of fluid site that sends, inverse direction].
                    streamingIndicesForReceivedDistributions[sharedSitesSeen] = contigSiteId
                                                                                * latticeInfo.GetNumVectors() + latticeInfo.GetInverseIndex(l);
                    ++sharedSitesSeen;
                }

            }

        }

    SiteRankIndex Domain::GetRankIndexFromGlobalCoords(const util::Vector3D<site_t> &globalSiteCoords) const {
        // Block identifiers (i, j, k) of the site (site_i, site_j, site_k)
        Vec16 blockCoords, localSiteCoords;
        GetBlockAndLocalSiteCoords(globalSiteCoords, blockCoords, localSiteCoords);
        return rank_for_site_store->GetSiteData(
                blockCoords,
                GetLocalSiteIdFromLocalSiteCoords(localSiteCoords)
        );
    }
        proc_t Domain::GetProcIdFromGlobalCoords(
                const util::Vector3D<site_t>& globalSiteCoords) const
        {
            return GetRankIndexFromGlobalCoords(globalSiteCoords)[0];
        }

        bool Domain::IsValidBlock(site_t i, site_t j, site_t k) const
        {
            if (i < 0 || i >= blockCounts.x())
                return false;
            if (j < 0 || j >= blockCounts.y())
                return false;
            if (k < 0 || k >= blockCounts.z())
                return false;

            return true;
        }

        bool Domain::IsValidBlock(const Vec16& blockCoords) const
        {
            for (int i =0 ; i < 3; ++i)
                if (blockCoords[i] < 0 || blockCoords[i] >= blockCounts[i])
                    return false;
            return true;
        }

        bool Domain::IsValidLatticeSite(const util::Vector3D<site_t>& siteCoords) const
        {
            using V = util::Vector3D<site_t>;
            return siteCoords.IsInRange(V::Zero(), sites - V::Ones());
        }

        site_t Domain::GetContiguousSiteId(util::Vector3D<site_t> location) const
        {
            // Block identifiers (i, j, k) of the site (site_i, site_j, site_k)
            Vec16 blockCoords, localSiteCoords;
            GetBlockAndLocalSiteCoords(location, blockCoords, localSiteCoords);
            // Pointer to the block
            const Block& lBlock = GetBlock(blockCoords);
            // Return pointer to site_data[site]
            return lBlock.GetLocalContiguousIndexForSite(GetLocalSiteIdFromLocalSiteCoords(localSiteCoords));
        }

        bool Domain::GetContiguousSiteId(const util::Vector3D<site_t>& globalLocation,
                                         proc_t& procId, site_t& siteId) const
        {
            if (!IsValidLatticeSite(globalLocation))
                return false;
            // convert global coordinates to local coordinates - i.e.
            // to location of block and location of site within block
            Vec16 blockCoords, localSiteCoords;
            GetBlockAndLocalSiteCoords(globalLocation, blockCoords, localSiteCoords);

            // get information for the block using the block location
            const Block& block = GetBlock(blockCoords);
            if (block.IsEmpty())
                return false;

            // get the local site id, i.e. its index within the block
            site_t localSiteIndex = GetLocalSiteIdFromLocalSiteCoords(localSiteCoords);

            // get the rank of the processor that owns the site
            procId = block.GetProcessorRankForSite(localSiteIndex);
            if (procId != comms.Rank())
                return false;
            if (procId == SITE_OR_BLOCK_SOLID) // means that the site is solid
                return false;

            // we only know enough information to determine solid/fluid for local sites
            // get the local contiguous index of the fluid site
            if (block.SiteIsSolid(localSiteIndex))
                return false;
            else
                siteId = block.GetLocalContiguousIndexForSite(localSiteIndex);

            return true;
        }

        util::Vector3D<site_t> Domain::GetGlobalCoords(
                site_t blockNumber, const util::Vector3D<site_t>& localSiteCoords) const
        {
            auto blockCoords = GetBlockIJK(blockNumber).as<site_t>();
            return GetGlobalCoords(blockCoords, localSiteCoords);
        }

        util::Vector3D<site_t> Domain::GetSiteCoordsFromSiteId(site_t siteId) const
        {
            util::Vector3D<site_t> siteCoords;
            siteCoords.z() = siteId % blockSize;
            site_t siteIJData = siteId / blockSize;
            siteCoords.y() = siteIJData % blockSize;
            siteCoords.x() = siteIJData / blockSize;
            return siteCoords;
        }

        void Domain::GetBlockAndLocalSiteCoords(const util::Vector3D<site_t>& location,
                                                Vec16& blockCoords,
                                                Vec16& siteCoords) const
        {
            blockCoords = (location / blockSize).as<U16>();
            siteCoords = (location % blockSize).as<U16>();
        }

    bool Domain::IsSiteDomainEdge(int rank, site_t local_idx) const {
        if (!remote_counts_cache.contains(rank)) {
            auto &tmp = remote_counts_cache[rank];
            shared_site_type_indices.Get(std::span{tmp}, rank);
        }
        auto const& counts = remote_counts_cache[rank];
        // Mid-domain indices are first
        auto n_mid_domain = counts[COLLISION_TYPES];
        return local_idx >= n_mid_domain;
    }

    Vec16 Domain::GetBlockIJK(site_t block) const
    {
        auto& t = rank_for_site_store->GetTree();
        return t.GetLeafCoords(block);
    }

        void Domain::Report(reporting::Dict& dictionary)
        {
            dictionary.SetIntValue("SITES", GetTotalFluidSites());
            dictionary.SetIntValue("BLOCKS", blockCount);
            dictionary.SetIntValue("SITESPERBLOCK", sitesPerBlockVolumeUnit);
            for (std::size_t n = 0; n < fluidSitesOnEachProcessor.size(); n++)
            {
                reporting::Dict proc = dictionary.AddSectionDictionary("PROCESSOR");
                proc.SetIntValue("RANK", n);
                proc.SetIntValue("SITES", fluidSitesOnEachProcessor[n]);
            }
        }
        neighbouring::NeighbouringDomain &Domain::GetNeighbouringData()
        {
            return *neighbouringData;
        }
        neighbouring::NeighbouringDomain const & Domain::GetNeighbouringData() const
        {
            return *neighbouringData;
        }

        int Domain::GetLocalRank() const
        {
            return comms.Rank();
        }

}
