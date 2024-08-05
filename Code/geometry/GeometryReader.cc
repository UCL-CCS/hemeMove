// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <cmath>
#include <list>
#include <algorithm>
#include <utility>
#include <zlib.h>

#include "io/formats/geometry.h"
#include "io/readers/XdrMemReader.h"
#include "geometry/decomposition/BasicDecomposition.h"
#include "geometry/decomposition/OptimisedDecomposition.h"
#include "geometry/GeometryReader.h"
#include "geometry/LookupTree.h"
#include "net/net.h"
#include "net/SparseExchange.h"
#include "net/IOCommunicator.h"
#include "log/Logger.h"
#include "util/span.h"
#include "util/numerical.h"
#include "util/Iterator.h"
#include "constants.h"

namespace hemelb::geometry
{
    namespace fmt = io::formats;
    using gmy = fmt::geometry;

    // Helper for checking that integers are allowed values for enums.
    template <typename Enum, Enum... allowed>
    struct EnumValidator {
    private:
        using INT = std::underlying_type_t<Enum>;

        // Want to turn the parameter pack into something iterable at
        // constexpr time, such as an initializer_list
        static constexpr bool IsValid(INT raw, std::initializer_list<Enum> vals) {
            for (auto& val: vals) {
                if (val == static_cast<Enum>(raw))
                    return true;
            }
            return false;
        }

    public:
        static Enum Run(INT raw) {
            if (IsValid(raw, {allowed...})) {
                return static_cast<Enum>(raw);
            } else {
                throw Exception() << "Invalid value for enum " << raw;
            }
        }
    };

    using SiteTypeValidator =  EnumValidator<gmy::SiteType,
            gmy::SiteType::SOLID,
            gmy::SiteType::FLUID>;
    using CutTypeValidator = EnumValidator<gmy::CutType,
            gmy::CutType::NONE,
            gmy::CutType::WALL,
            gmy::CutType::INLET,
            gmy::CutType::OUTLET>;
    using WallNormalAvailabilityValidator = EnumValidator<gmy::WallNormalAvailability,
            gmy::WallNormalAvailability::NOT_AVAILABLE,
            gmy::WallNormalAvailability::AVAILABLE>;

    GeometryReader::GeometryReader(const lb::LatticeInfo& latticeInfo,
                                   reporting::Timers &atimings, net::IOCommunicator ioComm) :
            latticeInfo(latticeInfo), computeComms(std::move(ioComm)), timings(atimings)
    {
    }

    GeometryReader::~GeometryReader()
    = default;

    GmyReadResult GeometryReader::LoadAndDecompose(const std::string& dataFilePath)
    {
        timings.fileRead().Start();

        // Open the file for read on node leaders
        if (computeComms.AmNodeLeader()) {
            file = net::MpiFile::Open(
                computeComms.GetLeadersComm(), dataFilePath, MPI_MODE_RDONLY, MPI_INFO_NULL
            );
        }

        log::Logger::Log<log::Debug, log::OnePerCore>("Reading file preamble");
        GmyReadResult geometry = ReadPreamble();

        log::Logger::Log<log::Debug, log::OnePerCore>("Reading file header");
        ReadHeader(geometry.GetBlockCount());
        timings.fileRead().Stop();

        {
            timings.initialDecomposition().Start();
            principalProcForEachBlock.resize(geometry.GetBlockCount());

            log::Logger::Log<log::Info, log::Singleton>("Creating block-level octree");
            auto blockTree = octree::build_block_tree(
                    geometry.GetBlockDimensions().as<octree::U16>(),
                    fluidSitesOnEachBlock
            );
            nFluidBlocks = blockTree.levels.back().node_ids.size();
            log::Logger::Log<log::Info, log::Singleton>(
                "Geometry has %lu / %ld active blocks, total %ld sites",
                nFluidBlocks,
                geometry.GetBlockCount(),
                blockTree.levels[0].sites_per_node[0]
            );

            // Get an initial base-level decomposition of the domain macro-blocks over processors.
            // This will later be improved upon by ParMetis.
            log::Logger::Log<log::Info, log::Singleton>("Beginning initial decomposition");
            decomposition::BasicDecomposition basicDecomposer(geometry,
                                                              computeComms.Size());

            // This vector only has entries for blocks that have a least one fluid site
            procForBlockOct = basicDecomposer.Decompose(blockTree, principalProcForEachBlock);
            geometry.block_store = std::make_unique<octree::DistributedStore>(
                    geometry.GetSitesPerBlock(),
                    std::move(blockTree),
                    procForBlockOct,
                    computeComms
            );
            if constexpr (build_info::VALIDATE_GEOMETRY) {
                log::Logger::Log<log::Info, log::Singleton>("Validating initial decomposition");
                basicDecomposer.Validate(principalProcForEachBlock, computeComms);
            }

            timings.initialDecomposition().Stop();
        }

        timings.fileRead().Start();
        {
          std::vector<U64> blocks_wanted;
          blocks_wanted.reserve((2*nFluidBlocks) / computeComms.Size());
          for (U64 i = 0; i < nFluidBlocks; ++i) {
            if (procForBlockOct[i] == computeComms.Rank())
              blocks_wanted.push_back(i);
          }
          ReadInBlocksWithHalo(geometry, blocks_wanted);
        }
        if constexpr (build_info::VALIDATE_GEOMETRY) {
            ValidateGeometry(geometry);
        }
        timings.fileRead().Stop();

        log::Logger::Log<log::Info, log::Singleton>("Optimising the domain decomposition.");
        timings.domainDecomposition().Start();

        // Having done an initial decomposition of the geometry, and read in the data, we optimise the
        // domain decomposition.
        log::Logger::Log<log::Debug, log::OnePerCore>("Beginning domain decomposition optimisation");
        OptimiseDomainDecomposition(geometry, principalProcForEachBlock);
        log::Logger::Log<log::Debug, log::OnePerCore>("Ending domain decomposition optimisation");

        if constexpr (build_info::VALIDATE_GEOMETRY) {
            log::Logger::Log<log::Info, log::Singleton>("Validating optimised decomposition");
            ValidateGeometry(geometry);
        }

        timings.domainDecomposition().Stop();

        return geometry;
    }

    std::vector<std::byte> GeometryReader::ReadAllProcesses(std::size_t start, unsigned nBytes)
    {
        // result
        std::vector<std::byte> buffer(nBytes);
        auto sp = to_span(buffer);
        if (computeComms.AmNodeLeader()) {
          // file is opened by the leaders comm only
          file.ReadAtAll(start, sp);
        }
        // Broadcast from node leader to others
        computeComms.GetNodeComm().Broadcast(sp, 0);
        return buffer;
    }

    /**
     * Read in the section at the beginning of the config file.
     */
    GmyReadResult GeometryReader::ReadPreamble()
    {
      auto preambleBuffer = ReadAllProcesses(0, gmy::PreambleLength);

      // Create an Xdr translator based on the read-in data.
      auto preambleReader = io::XdrMemReader(preambleBuffer.data(),
                                                           gmy::PreambleLength);

      uint32_t hlbMagicNumber, gmyMagicNumber, version;
      // Read in housekeeping values
      preambleReader.read(hlbMagicNumber);
      preambleReader.read(gmyMagicNumber);
      preambleReader.read(version);

      // Check the value of the HemeLB magic number.
      if (hlbMagicNumber != fmt::HemeLbMagicNumber)
      {
        throw Exception() << "This file does not start with the HemeLB magic number."
            << " Expected: " << unsigned(fmt::HemeLbMagicNumber)
            << " Actual: " << hlbMagicNumber;
      }

      // Check the value of the geometry file magic number.
      if (gmyMagicNumber != gmy::MagicNumber)
      {
        throw Exception() << "This file does not have the geometry magic number."
            << " Expected: " << unsigned(gmy::MagicNumber)
            << " Actual: " << gmyMagicNumber;
      }

      if (version != gmy::VersionNumber)
      {
        throw Exception() << "Version number incorrect."
            << " Supported: " << unsigned(gmy::VersionNumber)
            << " Input: " << version;
      }

      // Variables we'll read.
      // We use temporary vars here, as they must be the same size as the type in the file
      // regardless of the internal type used.
      uint32_t blocksX, blocksY, blocksZ, blockSize;

      auto read_check = [&](uint32_t& var) {
          preambleReader.read(var);
          if (var > std::uint32_t(std::numeric_limits<std::uint16_t>::max())) {
              throw Exception() << "size greater than 2^16 - 1!";
          }
      };
      // Read in the values.
      read_check(blocksX);
      read_check(blocksY);
      read_check(blocksZ);
      read_check(blockSize);

      // Read the padding unsigned int.
      unsigned paddingValue;
      preambleReader.read(paddingValue);

      return {Vec16(blocksX, blocksY, blocksZ), U16(blockSize)};
    }

    /**
     * Read the header section, with minimal information about each block.
     *
     * Results are placed in the member arrays fluidSitesPerBlock,
     * bytesPerCompressedBlock and bytesPerUncompressedBlock.
     */
    void GeometryReader::ReadHeader(site_t blockCount)
    {
      site_t headerByteCount = GetHeaderLength(blockCount);
      auto headerBuffer = ReadAllProcesses(gmy::PreambleLength, headerByteCount);

      // Create a Xdr translation object to translate from binary
      auto preambleReader = io::XdrMemReader(headerBuffer.data(),
                                             headerByteCount);

      fluidSitesOnEachBlock.reserve(blockCount);
      bytesPerCompressedBlock.reserve(blockCount);
      bytesPerUncompressedBlock.reserve(blockCount);

      // Read in all the data.
      for (site_t block = 0; block < blockCount; block++)
      {
        unsigned int sites, bytes, uncompressedBytes;
        preambleReader.read(sites);
        preambleReader.read(bytes);
        preambleReader.read(uncompressedBytes);

        fluidSitesOnEachBlock.push_back(sites);
        bytesPerCompressedBlock.push_back(bytes);
        bytesPerUncompressedBlock.push_back(uncompressedBytes);
      }
    }

    // Args are vectors with index representing a block ID.
    // First arg is which blocks this rank wants
    // Second arg is which rank "owns" them
    // Returns a pair with:
    // - a map from source rank to vector of block IDs wanted from it
    // - the total number of blocks wanted.
    auto ComputeSources(std::vector<bool> const& wanted, std::vector<int> const& ranks_they_are_on) {
      std::map<int, std::vector<std::size_t>> ans;
      auto const N = wanted.size();
      unsigned total = 0;
      for (std::size_t i = 0; i < N; ++i) {
        if (wanted[i]) {
          // This will default construct an empty vec if this is the first time
          auto& blocks = ans[ranks_they_are_on[i]];
          blocks.push_back(i);
          total += 1;
        }
      }
      return std::make_pair(std::move(ans), total);
    }


    // Return a map from block gmy index to compressed data, for
    // those blocks that where the predicate returns true
    template <std::predicate<std::size_t> PredT>
    auto GeometryReader::ReadCompressedBlockData(GmyReadResult const& geometry, PredT&& want) -> block_cache {
        // Strategy: read the entire file once and filter out those
        // blocks we want.

        // Going to read the file in large chunks (but note using
        // whole blocks) and only on the node leaders (using
        // collective MPI IO). Going to use the node communicator to
        // set up a shared memory allocation to hold this and then
        // filter in parallel.
        constexpr auto MiB = std::size_t(1) << 20;
        constexpr std::size_t MAX_GMY_BUFFER_SIZE = 64U * MiB;

        log::Logger::Log<log::Info, log::Singleton>("Streaming geometry data and caching required blocks.");
        log::Logger::Log<log::Debug, log::Singleton>("Maximum buffer size %lu B", MAX_GMY_BUFFER_SIZE);

        block_cache ans;

        // Work out where blocks live in the gmy **file**
        std::size_t const nBlocksGmy = bytesPerCompressedBlock.size();
        log::Logger::Log<log::Debug, log::Singleton>("Number of GMY blocks %lu", nBlocksGmy);
        std::size_t const dataStart = gmy::PreambleLength + GetHeaderLength(geometry.GetBlockCount());

        // N + 1 elements, elem i holds the start of block i, elem i+1 holds the end
        auto const blockBoundsGmy = [&] () {
            std::vector<std::size_t> ans(nBlocksGmy + 1);
            ans[0] = dataStart;
            std::inclusive_scan(
                bytesPerCompressedBlock.begin(), bytesPerCompressedBlock.end(),
                ans.begin() + 1, std::plus<std::size_t>(), dataStart
            );
            return ans;
        }();

        log::Logger::Log<log::Debug, log::Singleton>("Setup node level shared memory");

        MPI_Win win;
        std::byte* local_buf = nullptr;

        // Full size of buffer across node communicator
        auto total_buf_size = std::min(blockBoundsGmy[nBlocksGmy] - blockBoundsGmy[0],
                                       MAX_GMY_BUFFER_SIZE);

        auto&& nodeComm = computeComms.GetNodeComm();
        auto local_buf_size = (total_buf_size - 1) / nodeComm.Size() + 1;

        // Need contiguous memory.
        // Note local_buf has pointer into that process's "bit" of memory
        net::MpiCall{MPI_Win_allocate_shared}(
            local_buf_size, 1, MPI_INFO_NULL, nodeComm, &local_buf, &win
        );
        // Get the whole thing's start address
        std::byte* buf = local_buf - nodeComm.Rank() * local_buf_size;

        // Open a passive access epoch to the shared buffer
        net::MpiCall{MPI_Win_lock_all}(MPI_MODE_NOCHECK, win);

        // Get to work reading chunks
        std::size_t const* blockBoundsGmy_end = &*blockBoundsGmy.end();
        for (std::size_t i_first_block = 0; i_first_block < nBlocksGmy; /* end of loop */) {
          // Given we know which block we're starting at, figure out
          // the most whole blocks we can fit in the buffer.
          auto first_block_ptr = &blockBoundsGmy[i_first_block];
          auto max_read_pos = *first_block_ptr + total_buf_size;
          // upper bound gives the first elem after max_read_pos or _end if none
          auto end_ptr = std::upper_bound(first_block_ptr, blockBoundsGmy_end, max_read_pos) - 1;
          std::size_t n_blocks = end_ptr - first_block_ptr;
          auto read_size = *end_ptr - *first_block_ptr;
          log::Logger::Log<log::Debug, log::Singleton>("Reading blocks from %lu count %lu", i_first_block, n_blocks);

          // Only read on node leader
          if (computeComms.AmNodeLeader()) {
            auto sp = std::span<std::byte>(buf, read_size);
            // Collective on leaders comm
            file.ReadAtAll(*first_block_ptr, sp);
          }
          // Need to wait for leader to read
          nodeComm.Barrier();
          // Sync the memory
          net::MpiCall{MPI_Win_sync}(win);

          // Now we've read a chunk of the file. Go through it,
          // copying out the blocks we want.
          for (std::size_t i = 0; i < n_blocks; ++i) {
            auto block_gmy = i_first_block + i;
            if (want(block_gmy)) {
              // Recall std::map::operator[] will create the element.
              auto& block_data = ans[block_gmy];
              block_data.resize(bytesPerCompressedBlock[block_gmy]);

              auto buf_pos = blockBoundsGmy[block_gmy] - blockBoundsGmy[i_first_block];
              std::memcpy(block_data.data(), &buf[buf_pos], bytesPerCompressedBlock[block_gmy]);
            }
          }

          i_first_block += n_blocks;
        }
        // Close the access epoch
        net::MpiCall{MPI_Win_unlock_all}(win);
        // and free the window & buffer
        net::MpiCall{MPI_Win_free}(&win);

        log::Logger::Log<log::Debug, log::Singleton>("Finished caching blocks");
        return ans;
    }

    /**
     * Read and deserialise the necessary blocks from the file into
     * `geometry`. Collective.
     *
     * Initially, we have procForBlock which is the rank where each
     * block is wanted (or -1 if unknown).
     *
     * - Decide which blocks this rank want (i.e. procForBlock ==
     *   Rank() plus neigbours).
     *
     * - Collectively read the file, filtering blocks to ranks that
     *   want them.
     *
     * - Deserialise those blocks into the geometry.
     */
    void GeometryReader::ReadInBlocksWithHalo(GmyReadResult& geometry,
                                              const std::vector<U64>& blocksWanted)
    {
      // Create a list of which blocks to read in.
      timings.readBlocksPrelim().Start();

      // Populate the list of blocks to read (including a halo one block wide around all
      // local blocks).
      log::Logger::Log<log::Debug, log::OnePerCore>("Determining blocks to read");

      auto halo_wanted = DecideWhichBlocksToReadIncludingHalo(geometry, blocksWanted);

      // GMY file indexs of blocks we want.
      std::vector<std::size_t> wanted_gmys;
      wanted_gmys.reserve(halo_wanted.size());
      auto&& tree = geometry.block_store->GetTree();
      for (auto idx: halo_wanted) {
          auto ijk = tree.GetLeafCoords(idx);
          wanted_gmys.push_back(geometry.GetBlockIdFromBlockCoordinates(ijk));
      }
      // Recall we hit every block in GMY order, so sort this. We now
      // only have to test one value at a time and bump the iterator
      // forward when it matches.
      std::sort(wanted_gmys.begin(), wanted_gmys.end());

      auto compressed_block_data = ReadCompressedBlockData(
          geometry,
          [lower=wanted_gmys.cbegin(), upper=wanted_gmys.cend()] (std::size_t gmy) mutable {
            if (lower != upper && *lower == gmy) {
              ++lower;
              return true;
            } else {
              return false;
            }
          }
      );

      for (auto& [gmy_idx, data]: compressed_block_data) {
        DeserialiseBlock(geometry, data, gmy_idx);
      }
      timings.readBlocksPrelim().Stop();
    }

    void GeometryReader::DeserialiseBlock(
        GmyReadResult& geometry, std::vector<std::byte> const& compressedBlockData,
        site_t block_gmy
    ) {
        timings.readParse().Start();
        // Create an Xdr interpreter.
        auto blockData = DecompressBlockData(compressedBlockData,
                                             bytesPerUncompressedBlock[block_gmy]);
        io::XdrMemReader lReader(&blockData.front(), blockData.size());

        ParseBlock(geometry, block_gmy, lReader);

        // If debug-level logging, check that we've read in as many sites as anticipated.
        if constexpr (build_info::VALIDATE_GEOMETRY) {
          // Count the sites read,
          site_t numSitesRead = 0;
          for (site_t site = 0; site < geometry.GetSitesPerBlock(); ++site)
          {
            if (geometry.Blocks[block_gmy].Sites[site].targetProcessor != SITE_OR_BLOCK_SOLID)
            {
              ++numSitesRead;
            }
          }
          // Compare with the sites we expected to read.
          if (numSitesRead != fluidSitesOnEachBlock[block_gmy])
          {
            log::Logger::Log<log::Error, log::OnePerCore>("Was expecting %i fluid sites on block %i but actually read %i",
                                                          fluidSitesOnEachBlock[block_gmy],
                                                          block_gmy,
                                                          numSitesRead);
          }
        }
      timings.readParse().Stop();
    }

    std::vector<std::byte> GeometryReader::DecompressBlockData(const std::vector<std::byte>& compressed,
                                                          const unsigned int uncompressedBytes)
    {
      timings.unzip().Start();
      // For zlib return codes.
      int ret;

      // Set up the buffer for decompressed data. We know how long the the data is
      std::vector<std::byte> uncompressed(uncompressedBytes);

      // Set up the inflator
      z_stream stream;
      stream.zalloc = Z_NULL;
      stream.zfree = Z_NULL;
      stream.opaque = Z_NULL;
      stream.avail_in = compressed.size();
      stream.next_in = reinterpret_cast<unsigned char*>(const_cast<std::byte*>(compressed.data()));

      ret = inflateInit(&stream);
      if (ret != Z_OK)
        throw Exception() << "Decompression error for block";

      stream.avail_out = uncompressed.size();
      stream.next_out = reinterpret_cast<unsigned char*>(uncompressed.data());

      ret = inflate(&stream, Z_FINISH);
      if (ret != Z_STREAM_END)
        throw Exception() << "Decompression error for block";

      uncompressed.resize(uncompressed.size() - stream.avail_out);
      ret = inflateEnd(&stream);
      if (ret != Z_OK)
        throw Exception() << "Decompression error for block";

      timings.unzip().Stop();
      return uncompressed;
    }

    void GeometryReader::ParseBlock(GmyReadResult& geometry, const site_t block,
                                    io::XdrReader& reader)
    {
      // We start by clearing the sites on the block. We read the blocks twice (once before
      // optimisation and once after), so there can be sites on the block from the previous read.
      geometry.Blocks[block].Sites.clear();

      for (site_t localSiteIndex = 0; localSiteIndex < geometry.GetSitesPerBlock();
          ++localSiteIndex)
      {
        geometry.Blocks[block].Sites.push_back(ParseSite(reader));
      }
    }

    GeometrySite GeometryReader::ParseSite(io::XdrReader& reader)
    {
      // Read the site type
      unsigned readSiteType;
      bool success = reader.read(readSiteType);
      if (!success)
      {
        log::Logger::Log<log::Error, log::OnePerCore>("Error reading site type");
      }
      auto siteType = SiteTypeValidator::Run(readSiteType);
      GeometrySite readInSite(siteType == gmy::SiteType::FLUID);

      // If solid, there's nothing more to do.
      if (!readInSite.isFluid)
      {
        return readInSite;
      }

      // Prepare the links array to have enough space.
      readInSite.links.resize(latticeInfo.GetNumVectors() - 1);

      bool isGmyWallSite = false;

      // For each link direction...
      for (auto&& dir: gmy::Neighbourhood)
      {
        // read the type of the intersection and create a link...
        auto intersectionType = [&]() {
          unsigned readType;
          reader.read(readType);
          return CutTypeValidator::Run(readType);
        } ();

        GeometrySiteLink link;
        link.type = intersectionType;

        // walls have a floating-point distance to the wall...
        if (link.type == gmy::CutType::WALL)
        {
          isGmyWallSite = true;
          float distance;
          reader.read(distance);
          link.distanceToIntersection = distance;
        }
        // inlets and outlets (which together with none make up the other intersection types)
        // have an iolet id and a distance float...
        else if (link.type != gmy::CutType::NONE)
        {
          float distance;
          unsigned ioletId;
          reader.read(ioletId);
          reader.read(distance);

          link.ioletId = ioletId;
          link.distanceToIntersection = distance;
        }

        // Now, attempt to match the direction read from the local neighbourhood to one in the
        // lattice being used for simulation. If a match is found, assign the link to the read
        // site.
        for (Direction usedLatticeDirection = 1; usedLatticeDirection < latticeInfo.GetNumVectors();
            usedLatticeDirection++)
        {
          if (latticeInfo.GetVector(usedLatticeDirection) == dir)
          {
            // If this link direction is necessary to the lattice in use, keep the link data.
            readInSite.links[usedLatticeDirection - 1] = link;
            break;
          }
        }
      }

      auto normalAvailable = [&]() {
        unsigned normalAvailable;
        reader.read(normalAvailable);
        return WallNormalAvailabilityValidator::Run(normalAvailable);
      }();
      readInSite.wallNormalAvailable = (normalAvailable == gmy::WallNormalAvailability::AVAILABLE);

      if (readInSite.wallNormalAvailable != isGmyWallSite)
      {
        std::string msg = isGmyWallSite ?
          "wall fluid site without" :
          "bulk fluid site with";
        throw Exception() << "Malformed GMY file, " << msg
            << " a defined wall normal currently not allowed.";
      }

      if (readInSite.wallNormalAvailable)
      {
        reader.read(readInSite.wallNormal[0]);
        reader.read(readInSite.wallNormal[1]);
        reader.read(readInSite.wallNormal[2]);
      }

      return readInSite;
    }

    /**
     * This function is only called if in geometry-validation mode.
     * @param geometry
     */
    void GeometryReader::ValidateGeometry(const GmyReadResult& geometry)
    {
      log::Logger::Log<log::Debug, log::OnePerCore>("Validating the GlobalLatticeData");

      auto const SPB = geometry.GetSitesPerBlock();
      auto const NV = latticeInfo.GetNumVectors();
      constexpr auto UMAX = std::numeric_limits<unsigned>::max();

      std::vector<proc_t> myProcForSite(SPB);
      std::vector<unsigned> dummySiteData(SPB * NV);

      // We check the isFluid property and the link type for each direction
      // We also validate that each processor has the same beliefs about each site.
      for (site_t block_gmy = 0; block_gmy < geometry.GetBlockCount(); ++block_gmy) {
        auto const& block = geometry.Blocks[block_gmy];

        if (block.Sites.empty()) {
          std::fill(myProcForSite.begin(), myProcForSite.end(), SITE_OR_BLOCK_SOLID);
          std::fill(dummySiteData.begin(), dummySiteData.end(), UMAX);
        } else {
          for (site_t localSite = 0; localSite < SPB; ++localSite) {
            auto const& site = block.Sites[localSite];
            myProcForSite[localSite] = site.targetProcessor;

            auto dsd_start = localSite*NV;
            dummySiteData[dsd_start] = site.isFluid;
            for (Direction direction = 1; direction < NV; ++direction) {
              dummySiteData[dsd_start + direction] = site.isFluid ? unsigned(site.links[direction-1].type) : UMAX;
            }
          }
        }

        // Reduce using a maximum to find the actual processor for each site (ignoring the
        // invalid entries).
        std::vector<proc_t> procForSiteRecv = computeComms.AllReduce(myProcForSite, MPI_MAX);
        std::vector<unsigned> siteDataRecv = computeComms.AllReduce(dummySiteData, MPI_MIN);

        for (site_t site = 0; site < SPB; ++site) {
          if (myProcForSite[site] < 0)
            continue;

          if (myProcForSite[site] != procForSiteRecv[site]) {
            log::Logger::Log<log::Critical, log::OnePerCore>(
                "Site %li of block %li believed to be on %li but other process thinks %li",
                site, block_gmy, myProcForSite[site], procForSiteRecv[site]
            );
          }

          if (dummySiteData[site * NV] != siteDataRecv[site * NV]) {
            log::Logger::Log<log::Critical, log::OnePerCore>("Different fluid state was found for site %li on block %li. One: %li, Two: %li .",
                                                             site,
                                                             block_gmy,
                                                             dummySiteData[site*NV],
                                                             siteDataRecv[site*NV]
                                                             );
          }

          for (Direction dir = 1; dir < NV; ++dir) {
            if (dummySiteData[site * NV + dir] != siteDataRecv[site * NV+ dir]) {
              log::Logger::Log<log::Critical, log::OnePerCore>("Different link type was found for site %li, link %i on block %li. One: %li, Two: %li .",
                                                               site,
                                                               dir,
                                                               block_gmy,
                                                               dummySiteData[site*NV + dir],
                                                               siteDataRecv[site *NV + dir]);
            }
          }

        }
      }
    }


    // Go through blocks_wanted and add any 26-neighbouring blocks that are non-solid, using OCT ids.
    std::vector<U64> GeometryReader::DecideWhichBlocksToReadIncludingHalo(
        const GmyReadResult& geometry, const std::vector<U64>& blocks_wanted
    ) const {
      // Going to have to go from "compressed" octree order index
      // (idx) -> 3D grid coordinate (ijk) to compute neighbours. Get
      // set up for this.
      auto const block_dims = geometry.GetBlockDimensions();
      constexpr auto U16_MAX = std::numeric_limits<U16>::max();
      auto&& tree = geometry.block_store->GetTree();

      // Start with no blocks wanted.
      std::vector<bool> want_block_on_this_rank(nFluidBlocks, false);

      // Main loop
      for (auto block_idx: blocks_wanted) {
          // Could set this here, but will cover in loop over halo below
          // want_block_on_this_rank[block_idx] = true;

          auto block_ijk = tree.GetLeafCoords(block_idx);

          // Compute the bounds of neighbours to consider
          Vec16 lo, hi;
          for (int d = 0; d < 3; ++d) {
              // Beware the "usual arithmetic conversions"!
              U16 const v = block_ijk[d];
              lo[d] = v > 0U ? v - 1U : 0U;
              // Note, we will use inclusive hi limit below so max
              // allowed value is block_dims[d] - 1
              U16 const w = v < U16_MAX ? v + 1U : v;
              hi[d] = w < block_dims[d] ? w : block_dims[d] - 1;
          }

          // 3D loop
          for (U16 ni = lo[0]; ni <= hi[0]; ++ni)
              for (U16 nj = lo[1]; nj <= hi[1]; ++nj)
                  for (U16 nk = lo[2]; nk <= hi[2]; ++nk) {
                      auto neigh_ijk = Vec16{ni, nj, nk};
                      auto neigh_idx = tree.GetPath(neigh_ijk).leaf();
                      // Recall that the octree will return a path
                      // with "no child" for all levels where the
                      // requested node doesn't exist.
                      if (neigh_idx != octree::Level::NC)
                          want_block_on_this_rank[neigh_idx] = true;
                  }
      }

      std::vector<U64> ans;
      for (unsigned i = 0; i < nFluidBlocks; ++i) {
        if (want_block_on_this_rank[i])
          ans.push_back(i);
      }
      return ans;
    }

    void GeometryReader::OptimiseDomainDecomposition(GmyReadResult& geometry,
                                                     const std::vector<proc_t>& procForEachBlock)
    {
      decomposition::OptimisedDecomposition optimiser(timings,
                                                      computeComms,
                                                      geometry,
                                                      latticeInfo);

      timings.reRead().Start();
      log::Logger::Log<log::Debug, log::OnePerCore>("Rereading blocks");
      // Reread the blocks based on the ParMetis decomposition.
      RereadBlocks(geometry,
                   optimiser.GetStaying(),
                   optimiser.GetArriving());
      timings.reRead().Stop();

      timings.moves().Start();
      // Implement the decomposition now that we have read the necessary data.
      log::Logger::Log<log::Debug, log::OnePerCore>("Implementing moves");
      ImplementMoves(geometry,
                     optimiser.GetStaying(),
                     optimiser.GetArriving(),
                     optimiser.GetLeaving());
      timings.moves().Stop();
    }

    // The header section of the config file contains a number of records.
    site_t GeometryReader::GetHeaderLength(site_t blockCount) const
    {
      return gmy::HeaderRecordLength * blockCount;
    }

    // Iterator to advance through the moves vector to the first move
    // in the next block, recalling that moves are sorted by block and
    // there's a max number of sites per block.
    struct only_block_id_iterator {
        using Iter = SiteVec::const_iterator;

        site_t max_step;
        Iter pos;
        Iter end;

        // Only care about block ID
        auto operator*() const {
          return (*pos)[0];
        }

        only_block_id_iterator& operator++() {
            auto block = **this;
            ++pos; // Guard the case where have a fully fluid block
            auto limit = std::min(pos + max_step, end);
            pos = std::upper_bound(
                pos, limit,
                block,
                [](U64 l, SiteDesc const& r) {return l < r[0]; }
            );
            return *this;
        }

        operator bool() const {
            return pos < end;
        }
    };

    void GeometryReader::RereadBlocks(GmyReadResult& geometry, SiteVec const& staying,
                                      MovesMap const& arriving)
    {
      // Go through the sites that we will end up with, and compute
      // the unique list of blocks (by OCT index) that we need. Recall
      // that the moves lists are sorted by block and then site id
      // (tho we only care about block).
      std::vector<U64> optimal_blocks;

      // We can be simple for the staying blocks as opt is empty at
      // start and know:
      // - that the blocks monotonically increase thru the array
      // - that there are at most SPB sites per block
      auto const SPB = geometry.GetSitesPerBlock();
      for (auto it = only_block_id_iterator{SPB, staying.begin(), staying.end()}; it; ++it) {
          optimal_blocks.push_back(*it);
      }

      // For the arriving sites, need to deal with insertion
      for (auto& [src_rank, data]: arriving) {
        // Recall that each rank's data is sorted by block ID, thus we
        // can search for the insert point more efficiently.
        auto opt_pos = optimal_blocks.begin();

        for (auto it = only_block_id_iterator{SPB, data.begin(), data.end()}; it; ++it) {
          auto block = *it;

          auto lb = std::lower_bound(opt_pos, optimal_blocks.end(), block);
          // BUT adding a block may invalidate the opt_pos iterator
          if (lb ==  optimal_blocks.end()) {
            // Not found and greater than any value
            optimal_blocks.push_back(block);
            opt_pos = optimal_blocks.end();
          } else if (*lb == block) {
            // already there
            ++opt_pos;
          } else {
            // Not found and need to insert and move past the inserted one
            opt_pos = ++optimal_blocks.insert(lb, block);
          }
        }
      }

      // Reread the blocks into the GlobalLatticeData now.
      ReadInBlocksWithHalo(geometry, optimal_blocks);
    }

    void GeometryReader::ImplementMoves(GmyReadResult& geometry,
                                        SiteVec const& staying,
                                        MovesMap const& arriving,
                                        MovesMap const& leaving) const
    {
        // Given a vector of sites (sorted by block, then site index
        // within the block), assign the process to that site in the
        // GmyReadResult.
        auto set_rank_for_sites = [&] (SiteVec const& sites, int rank) {

            auto block_start = sites.begin();
            auto end = sites.end();
            auto&& tree = geometry.block_store->GetTree();

            while (block_start != end) {
                auto const block_oct = (*block_start)[0];
                site_t const nsites = tree.levels[tree.n_levels].sites_per_node[block_oct];
                auto const block_end = (++only_block_id_iterator{nsites, block_start, end}).pos;

                auto const block_ijk = tree.GetLeafCoords(block_oct);
                auto const block_gmy = geometry.GetBlockIdFromBlockCoordinates(block_ijk);

                auto& block_sites = geometry.Blocks[block_gmy].Sites;
                for (auto it = block_start; it < block_end; ++it) {
                    auto [block, site_idx] = *it;
                    HASSERT(block == block_oct);
                    HASSERT(!block_sites.empty());

                    auto& site = block_sites[site_idx];
                    HASSERT(site.isFluid);
                    site.targetProcessor = rank;
                }

                block_start = block_end;
            }
        };

        // First set proc for staying sites
        set_rank_for_sites(staying, computeComms.Rank());
        // Then arriving sites
        for (auto const& [_, sites]: arriving) {
          set_rank_for_sites(sites, computeComms.Rank());
        }
        // We could set the leaving sites, but can't guarantee we get
        // all the neighbours of our sites, so no saving of time
        // really.
    }

}
