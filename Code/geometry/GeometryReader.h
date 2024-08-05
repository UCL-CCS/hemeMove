// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_GEOMETRY_GEOMETRYREADER_H
#define HEMELB_GEOMETRY_GEOMETRYREADER_H

#include <vector>
#include <string>

#include "io/readers/XdrReader.h"
#include "lb/lattices/LatticeInfo.h"
#include "lb/LbmParameters.h"
#include "net/mpi.h"
#include "geometry/ParmetisForward.h"
#include "reporting/Timers.h"
#include "util/Vector3D.h"
#include "units.h"
#include "geometry/GmyReadResult.h"
#include "geometry/needs/Needs.h"

#include "net/MpiFile.h"

namespace hemelb::geometry
{

    class GeometryReader
    {
    public:
        using BlockLocation = util::Vector3D<site_t>;

        GeometryReader(const lb::LatticeInfo&,
                       reporting::Timers &timings, net::IOCommunicator ioComm);
        ~GeometryReader();

        GmyReadResult LoadAndDecompose(const std::string& dataFilePath);

    private:
        // Read from the file into a buffer on all processes.
        // This is collective and start and nBytes must be the same on all ranks.
        std::vector<std::byte> ReadAllProcesses(std::size_t startBytes, unsigned nBytes);

        // Read the preamble and create the empty read result.
        GmyReadResult ReadPreamble();

        // Read the block header with basic size data for each block in
        // the domain.
        void ReadHeader(site_t blockCount);

        // map from block GMY index to vect of data
        using block_cache = std::map<U64, std::vector<std::byte>>;

        // Load compressed data from the file if the predicate is true for that block's GMY index.
        template <std::predicate<std::size_t> PredT>
        block_cache ReadCompressedBlockData(GmyReadResult const& geometry, PredT&& p);


        // Given a vector of the block OCT ids that we want, add a
        // one-block halo and read those blocks into the
        // GmyReadResult.
        void ReadInBlocksWithHalo(GmyReadResult& geometry,
                                  const std::vector<U64>& blocksWanted);

        // Given a vector of the block OCT ids wanted, add a one-block
        // halo and return a vector of those block OCT ids.
        std::vector<U64> DecideWhichBlocksToReadIncludingHalo(
            const GmyReadResult& geometry, const std::vector<U64>& blocks_wanted
        ) const;

        // Parse a compressed block of data into the geometry at the given GMY index
        void DeserialiseBlock(GmyReadResult& geometry,
                              std::vector<std::byte> const& compressed_data, site_t block_gmy);

        // Decompress the block data
        std::vector<std::byte> DecompressBlockData(const std::vector<std::byte>& compressed,
                                              const unsigned int uncompressedBytes);

        // Given a reader for a block's data, parse that into the
        // GmyReadResult at the given index.
        void ParseBlock(GmyReadResult& geometry, const site_t block,
                        io::XdrReader& reader);

        // Parse the next site from the XDR reader
        GeometrySite ParseSite(io::XdrReader& reader);

        // Use the OptimisedDecomposition class to refine a simple,
        // block-level initial decomposition.
        void OptimiseDomainDecomposition(GmyReadResult& geometry,
                                         const std::vector<proc_t>& procForEachBlock);

        // Check for self-consistency
        void ValidateGeometry(const GmyReadResult& geometry);

        // Get the length of the header section, given the number of blocks.
        site_t GetHeaderLength(site_t blockCount) const;

        // Given knowledge of which sites are needed, read the blocks
        // required (i.e. the blocks the sites belong to plus a
        // one-block halo) into the GmyReadResult.
        void RereadBlocks(GmyReadResult& geometry, SiteVec const& staying_sites,
                          MovesMap const& arriving_sites);

        // Set the rank of sites assigned to this process in the GmyReadResult.
        void ImplementMoves(GmyReadResult& geometry, SiteVec const& staying,
                            MovesMap const& arriving, MovesMap const& leaving) const;

        //! Info about the connectivity of the lattice.
        const lb::LatticeInfo& latticeInfo;

        // File accessed to read in the geometry data.
        //
        // We are going
        // to read the whole file collectively and sequentially in
        // chunks. The node leaders will read data either into private
        // memory and broadcast or into node-shared memory and do a
        // barrier/fence.
        net::MpiFile file;

        //! Communicator for all ranks that will need a slice of the geometry
        net::IOCommunicator computeComms;

        //! How many blocks with at least one fluid site
        U64 nFluidBlocks;
        //! The number of fluid sites on each block in the file.
        std::vector<site_t> fluidSitesOnEachBlock;
        //! The number of bytes each block in the file takes up while still compressed.
        std::vector<unsigned int> bytesPerCompressedBlock;
        //! The number of bytes each block in the file takes up when uncompressed.
        std::vector<unsigned int> bytesPerUncompressedBlock;
        //! The process assigned to each block.
        std::vector<proc_t> principalProcForEachBlock;
        //! The process for fluid-containing blocks in octree order
        std::vector<proc_t> procForBlockOct;

        //! Timings object for recording the time taken for each step of the domain decomposition.
        hemelb::reporting::Timers &timings;
    };
}

#endif /* HEMELB_GEOMETRY_GEOMETRYREADER_H */
