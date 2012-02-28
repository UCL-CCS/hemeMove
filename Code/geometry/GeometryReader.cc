#include <cmath>
#include <list>
#include <map>
#include <zlib.h>

#include "debug/Debugger.h"
#include "io/formats/geometry.h"
#include "io/writers/xdr/XdrMemReader.h"
#include "geometry/GeometryReader.h"
#include "lb/lattices/D3Q27.h"
#include "net/net.h"
#include "topology/NetworkTopology.h"
#include "log/Logger.h"
#include "util/utilityFunctions.h"
#include "constants.h"

namespace hemelb
{
  namespace geometry
  {

    GeometryReader::GeometryReader(const bool reserveSteeringCore, GeometryReadResult& readResult, reporting::Timers &atimings) :
        readingResult(readResult), timings(atimings)
    {
      // Get the group of all procs.
      MPI_Group worldGroup;
      MPI_Comm_group(MPI_COMM_WORLD, &worldGroup);

      participateInTopology = !reserveSteeringCore || topology::NetworkTopology::Instance()->GetLocalRank() != 0;

      // Create our own group, without the root node.
      if (reserveSteeringCore && topology::NetworkTopology::Instance()->GetProcessorCount() > 1)
      {
        int lExclusions[1] = { 0 };
        MPI_Group_excl(worldGroup, 1, lExclusions, &topologyGroup);
      }
      else
      {
        topologyGroup = worldGroup;
      }

      // Create a communicator just for the domain decomposition.
      MPI_Comm_create(MPI_COMM_WORLD, topologyGroup, &topologyComm);

      // Each rank needs to know its rank wrt the domain
      // decomposition.
      if (participateInTopology)
      {
        int temp = 0;
        MPI_Comm_rank(topologyComm, &topologyRank);
        MPI_Comm_size(topologyComm, &temp);
        topologySize = (unsigned int) temp;
      }
      else
      {
        topologyRank = -1;
        topologySize = 0;
      }
    }

    GeometryReader::~GeometryReader()
    {
      MPI_Group_free(&topologyGroup);

      // Note that on rank 0, this is the same as MPI_COMM_WORLD.
      if (participateInTopology)
      {
        MPI_Comm_free(&topologyComm);
      }
    }

    void GeometryReader::LoadAndDecompose(std::string& dataFilePath)
    {
      timings[hemelb::reporting::Timers::fileRead].Start();

      int error;

      // Open the file using the MPI parallel I/O interface at the path
      // given, in read-only mode.
      MPI_Info_create(&fileInfo);

      // Create hints about how we'll read the file. See Chapter 13, page 400 of the MPI 2.2 spec.
      std::string accessStyle = "access_style";
      std::string accessStyleValue = "sequential";
      std::string buffering = "collective_buffering";
      std::string bufferingValue = "true";

      MPI_Info_set(fileInfo, &accessStyle[0], &accessStyleValue[0]);
      MPI_Info_set(fileInfo, &buffering[0], &bufferingValue[0]);

      // Open the file.
      error = MPI_File_open(MPI_COMM_WORLD, &dataFilePath[0], MPI_MODE_RDONLY, fileInfo, &file);

      currentCommRank = topology::NetworkTopology::Instance()->GetLocalRank();
      currentCommSize = topology::NetworkTopology::Instance()->GetProcessorCount();
      currentComm = MPI_COMM_WORLD;

      if (error != 0)
      {
        log::Logger::Log<log::Info, log::OnePerCore>("Unable to open file %s, exiting", dataFilePath.c_str());
        fflush(0x0);
        exit(0x0);
      }
      else
      {
        log::Logger::Log<log::Debug, log::OnePerCore>("Opened config file %s", dataFilePath.c_str());
      }
      fflush(NULL);

      // Set the view to the file.
      std::string mode = "native";
      MPI_File_set_view(file, 0, MPI_CHAR, MPI_CHAR, &mode[0], fileInfo);

      log::Logger::Log<log::Debug, log::OnePerCore>("Reading file preamble");
      ReadPreamble();

      // Read the file header.
      log::Logger::Log<log::Debug, log::OnePerCore>("Reading file header");

      fluidSitesPerBlock = new site_t[readingResult.GetBlockCount()];
      bytesPerCompressedBlock = new unsigned int[readingResult.GetBlockCount()];
      bytesPerUncompressedBlock = new unsigned int[readingResult.GetBlockCount()];
      procForEachBlock = new proc_t[readingResult.GetBlockCount()];

      ReadHeader();
      timings[hemelb::reporting::Timers::initialDecomposition].Start();
      // Perform an initial decomposition, of which processor should read each block.
      log::Logger::Log<log::Debug, log::OnePerCore>("Beginning initial decomposition");

      if (!participateInTopology)
      {
        for (site_t block = 0; block < readingResult.GetBlockCount(); ++block)
        {
          procForEachBlock[block] = -1;
        }
      }
      else
      {
        BlockDecomposition();

        ValidateProcForEachBlock();
      }
      timings[hemelb::reporting::Timers::initialDecomposition].Stop();
      // Perform the initial read-in.
      log::Logger::Log<log::Debug, log::OnePerCore>("Reading in my blocks");

      // Close the file - only the ranks participating in the topology need to read it again.
      MPI_File_close(&file);

      if (participateInTopology)
      {
        // Reopen in the file just between the nodes in the topology decomposition. Read in blocks
        // local to this node.
        MPI_File_open(topologyComm, &dataFilePath[0], MPI_MODE_RDONLY, fileInfo, &file);

        currentCommRank = topologyRank;
        currentCommSize = topologySize;
        currentComm = topologyComm;

        ReadInLocalBlocks(procForEachBlock, topologyRank);

        if (log::Logger::ShouldDisplay<log::Debug>())
        {
          ValidateAllReadData();
        }
      }

      timings[hemelb::reporting::Timers::fileRead].Stop();

      hemelb::log::Logger::Log<hemelb::log::Warning, hemelb::log::Singleton>("Begin optimising the domain decomposition.");
      timings[hemelb::reporting::Timers::domainDecomposition].Start();
      // Optimise.
      if (participateInTopology)
      {
        log::Logger::Log<log::Debug, log::OnePerCore>("Beginning domain decomposition optimisation");
        OptimiseDomainDecomposition(procForEachBlock);
        log::Logger::Log<log::Debug, log::OnePerCore>("Ending domain decomposition optimisation");

        if (log::Logger::ShouldDisplay<log::Debug>())
        {
          ValidateAllReadData();
        }

        MPI_File_close(&file);
      }

      // Finish up - close the file, set the timings, deallocate memory.
      MPI_Info_free(&fileInfo);

      timings[hemelb::reporting::Timers::domainDecomposition].Stop();

      delete[] fluidSitesPerBlock;
      delete[] bytesPerCompressedBlock;
      delete[] bytesPerUncompressedBlock;
      delete[] procForEachBlock;
    }

    /**
     * Read in the section at the beginning of the config file.
     */
    void GeometryReader::ReadPreamble()
    {
      const unsigned preambleBytes = io::formats::geometry::PreambleLength;

      // Read in the file preamble into a buffer.
      char preambleBuffer[preambleBytes];

      if (currentCommRank == HEADER_READING_RANK)
      {
        MPI_File_read(file, preambleBuffer, preambleBytes, MpiDataType(preambleBuffer[0]), MPI_STATUS_IGNORE);
      }

      MPI_Bcast(preambleBuffer, preambleBytes, MpiDataType<char>(), HEADER_READING_RANK, currentComm);

      // Create an Xdr translator based on the read-in data.
      hemelb::io::writers::xdr::XdrReader preambleReader = hemelb::io::writers::xdr::XdrMemReader(preambleBuffer,
                                                                                                  preambleBytes);

      unsigned hlbMagicNumber, gmyMagicNumber, version;
      // Read in housekeeping values
      preambleReader.readUnsignedInt(hlbMagicNumber);
      preambleReader.readUnsignedInt(gmyMagicNumber);
      preambleReader.readUnsignedInt(version);

      // Check the value of the HemeLB magic number.
      if (hlbMagicNumber != io::formats::HemeLbMagicNumber)
      {
        log::Logger::Log<log::Info, log::OnePerCore>("This file starts with %d, not the HemeLB magic number %d.",
                                                     hlbMagicNumber,
                                                     io::formats::HemeLbMagicNumber);
        exit(1);
      }

      // Check the value of the geometry file magic number.
      if (gmyMagicNumber != io::formats::geometry::MagicNumber)
      {
        log::Logger::Log<log::Info, log::OnePerCore>("This file is not a geometry file: had %d, not the geometry magic number %d.",
                                                     gmyMagicNumber,
                                                     io::formats::geometry::MagicNumber);
        exit(1);
      }

      log::Logger::Log<log::Warning, log::OnePerCore>("Geometry file version: %d", version);

      // Variables we'll read.
      // We use temporary vars here, as they must be the same size as the type in the file
      // regardless of the internal type used.
      unsigned int blocksX, blocksY, blocksZ, blockSize;

      // Read in the values.
      preambleReader.readUnsignedInt(blocksX);
      preambleReader.readUnsignedInt(blocksY);
      preambleReader.readUnsignedInt(blocksZ);
      preambleReader.readUnsignedInt(blockSize);
      preambleReader.readDouble(readingResult.voxelSize);
      for (unsigned int i = 0; i < 3; ++i)
      {
        preambleReader.readDouble(readingResult.origin[i]);
      }

      // Read the padding unsigned int.
      unsigned paddingValue;
      preambleReader.readUnsignedInt(paddingValue);

      readingResult.blocks.x = blocksX;
      readingResult.blocks.y = blocksY;
      readingResult.blocks.z = blocksZ;

      readingResult.blockSize = blockSize;

      // Ensure the reading result container has enough space to contain all the blocks.
      readingResult.Blocks.resize(readingResult.GetBlockCount());
    }

    /**
     * Read the header section, with minimal information about each block.
     *
     * Results are placed in the member arrays fluidSitesPerBlock,
     * bytesPerCompressedBlock and bytesPerUncompressedBlock.
     */
    void GeometryReader::ReadHeader()
    {
      site_t headerByteCount = GetHeaderLength(readingResult.GetBlockCount());
      // Allocate a buffer to read into, then do the reading.
      char* headerBuffer = new char[headerByteCount];

      if (currentCommRank == HEADER_READING_RANK)
      {
        MPI_File_read(file, headerBuffer, (int) headerByteCount, MpiDataType(headerBuffer[0]), MPI_STATUS_IGNORE);
      }

      MPI_Bcast(headerBuffer, (int) headerByteCount, MpiDataType<char>(), HEADER_READING_RANK, currentComm);

      // Create a Xdr translation object to translate from binary
      hemelb::io::writers::xdr::XdrReader preambleReader =
          hemelb::io::writers::xdr::XdrMemReader(headerBuffer, (unsigned int) headerByteCount);

      // Read in all the data.
      for (site_t block = 0; block < readingResult.GetBlockCount(); block++)
      {
        unsigned int sites, bytes, uncompressedBytes;
        preambleReader.readUnsignedInt(sites);
        preambleReader.readUnsignedInt(bytes);
        preambleReader.readUnsignedInt(uncompressedBytes);

        fluidSitesPerBlock[block] = sites;
        bytesPerCompressedBlock[block] = bytes;
        bytesPerUncompressedBlock[block] = uncompressedBytes;
      }

      delete[] headerBuffer;
    }

    /**
     * Read in the necessary blocks from the file.
     * TODO: explain what this method does
     */
    void GeometryReader::ReadInLocalBlocks(const proc_t* unitForEachBlock,
                                           const proc_t localRank)
    {
      // Create a list of which blocks to read in.
      timings[hemelb::reporting::Timers::readBlocksPrelim].Start();
      bool* readBlock = new bool[readingResult.GetBlockCount()];
      log::Logger::Log<log::Debug, log::OnePerCore>("Determining blocks to read");
      DecideWhichBlocksToRead(readBlock, unitForEachBlock, localRank);

      if (log::Logger::ShouldDisplay<log::Debug>())
      {
        log::Logger::Log<log::Debug, log::OnePerCore>("Validating block sizes");
        for (site_t block = 0; block < readingResult.GetBlockCount(); ++block)
        {
          if (bytesPerUncompressedBlock[block] >
              io::formats::geometry::GetMaxBlockRecordLength(readingResult.blockSize,
                                                             fluidSitesPerBlock[block]))
          {
            log::Logger::Log<log::Debug, log::OnePerCore>("Block %i is %i bytes when the longest possible block should be %i bytes",
                                                          block,
                                                          bytesPerUncompressedBlock[block],
                                                          io::formats::geometry::GetMaxBlockRecordLength(readingResult.blockSize,
                                                                                                         fluidSitesPerBlock[block]));
          }
        }
      }

      log::Logger::Log<log::Debug, log::OnePerCore>("Informing reading cores of block needs");
      net::Net net = net::Net(currentComm);
      Decomposition *decomposition=new Decomposition(readingResult.GetBlockCount(),
        readBlock,util::NumericalFunctions::min(READING_GROUP_SIZE, currentCommSize),net,currentComm,currentCommRank,currentCommSize);

      
      timings[hemelb::reporting::Timers::readBlocksPrelim].Stop();
      log::Logger::Log<log::Debug, log::OnePerCore>("Reading blocks");
      timings[hemelb::reporting::Timers::readBlocksAll].Start();
      
       // Set the view and read in.
        MPI_Offset offset = io::formats::geometry::PreambleLength + GetHeaderLength(readingResult.GetBlockCount());
      // Track the next block we should look at.
      for (site_t nextBlockToRead = 0; nextBlockToRead < readingResult.GetBlockCount(); ++nextBlockToRead)
      {
        ReadInBlock(offset,
                    decomposition->ProcessorsNeedingBlock(nextBlockToRead),
                    nextBlockToRead,
                    fluidSitesPerBlock[nextBlockToRead],
                    bytesPerCompressedBlock[nextBlockToRead],
                    bytesPerUncompressedBlock[nextBlockToRead],
                    readBlock[nextBlockToRead]);

        offset += bytesPerCompressedBlock[nextBlockToRead];
      }

      delete decomposition;
      delete[] readBlock;
      timings[hemelb::reporting::Timers::readBlocksAll].Stop();
    }

    void GeometryReader::ReadInBlock(MPI_Offset offsetSoFar,
                                     const std::vector<proc_t>& procsWantingThisBlock,
                                     const site_t blockNumber,
                                     const site_t sites,
                                     const unsigned int compressedBytes,
                                     const unsigned int uncompressedBytes,
                                     const int neededOnThisRank)
    {
      // Easy case if there are no sites on the block.
      if (sites <= 0)
      {
        return;
      }
      std::vector<char> compressedBlockData;
      proc_t readingCore = GetReadingCoreForBlock(blockNumber);

      net::Net net = net::Net(currentComm);
      
      if (readingCore == currentCommRank)
      {
        timings[hemelb::reporting::Timers::readBlock].Start();
        // Read the data.
        compressedBlockData.resize(compressedBytes);
        MPI_File_read_at(file, offsetSoFar, &compressedBlockData.front(), compressedBytes, MPI_CHAR, MPI_STATUS_IGNORE);

        // Spread it.
        for (std::vector<proc_t>::const_iterator receiver = procsWantingThisBlock.begin(); receiver != procsWantingThisBlock.end(); receiver++)
        {
          if (*receiver != currentCommRank)
          {
            net.RequestSend(&compressedBlockData.front(), compressedBytes, *receiver);
          }
        }
        timings[hemelb::reporting::Timers::readBlock].Stop();
      }
      else if (neededOnThisRank)
      {
        compressedBlockData.resize(compressedBytes);
        net.RequestReceive(&compressedBlockData.front(), compressedBytes, readingCore);
      }
      else
      {
        return;
      }
      timings[hemelb::reporting::Timers::readNet].Start();
      net.Send();
      net.Receive();
      net.Wait();
      timings[hemelb::reporting::Timers::readNet].Stop();
      timings[hemelb::reporting::Timers::readParse].Start();
      if (neededOnThisRank)
      {
        // Create an Xdr interpreter.
        std::vector<char> blockData = DecompressBlockData(compressedBlockData, uncompressedBytes);
        io::writers::xdr::XdrMemReader lReader(&blockData.front(), blockData.size());

        ParseBlock(blockNumber, lReader);

        // If debug-level logging, check that we've read in as many sites as anticipated.
        if (log::Logger::ShouldDisplay<log::Debug>())
        {
          // Count the sites read,
          site_t numSitesRead = 0;
          for (site_t site = 0; site < readingResult.GetSitesPerBlock(); ++site)
          {
            if (readingResult.Blocks[blockNumber].Sites[site].targetProcessor != BIG_NUMBER2)
            {
              ++numSitesRead;
            }
          }
          // Compare with the sites we expected to read.
          if (numSitesRead != sites)
          {
            log::Logger::Log<log::Debug, log::OnePerCore>("Was expecting %i fluid sites on block %i but actually read %i",
                                                          sites,
                                                          blockNumber,
                                                          numSitesRead);
          }
        }
      }
      timings[hemelb::reporting::Timers::readParse].Stop();
    }

    std::vector<char> GeometryReader::DecompressBlockData(const std::vector<char>& compressed,
                                                          const unsigned int uncompressedBytes) {
      timings[hemelb::reporting::Timers::unzip].Start();                        
      // For zlib return codes.
      int ret;

      // Set up the buffer for decompressed data. We know how long the the data is
      std::vector<char> uncompressed(uncompressedBytes);

      // Set up the inflator
      z_stream stream;
      stream.zalloc = Z_NULL;
      stream.zfree = Z_NULL;
      stream.opaque = Z_NULL;
      stream.avail_in = compressed.size();
      stream.next_in = reinterpret_cast<unsigned char*>(const_cast<char*>(&compressed.front()));

      ret = inflateInit(&stream);
      if (ret != Z_OK) {
        log::Logger::Log<log::Debug, log::OnePerCore>("Decompression error for block");
        std::exit(1);
      }
      stream.avail_out = uncompressed.size();
      stream.next_out = reinterpret_cast<unsigned char*>(&uncompressed.front());

      ret = inflate(&stream, Z_FINISH);
      if (ret != Z_STREAM_END) {
        log::Logger::Log<log::Debug, log::OnePerCore>("Decompression error for block");
        std::exit(1);
      }

      uncompressed.resize(uncompressed.size() - stream.avail_out);
      ret = inflateEnd(&stream);
      if (ret != Z_OK) {
        log::Logger::Log<log::Debug, log::OnePerCore>("Decompression error for block");
        std::exit(1);
      }
      timings[hemelb::reporting::Timers::unzip].Stop();        
      return uncompressed;
    }

    void GeometryReader::ParseBlock(const site_t block, io::writers::xdr::XdrReader& reader)
    {
      readingResult.Blocks[block].Sites.clear();

      for (site_t localSiteIndex = 0; localSiteIndex < readingResult.GetSitesPerBlock(); ++localSiteIndex)
      {
        readingResult.Blocks[block].Sites.push_back(ParseSite(reader));
      }
    }

    SiteReadResult GeometryReader::ParseSite(io::writers::xdr::XdrReader& reader)
    {
      // Read the fluid property.
      unsigned isFluid;

      if (!reader.readUnsignedInt(isFluid))
      {
        log::Logger::Log<log::Info, log::OnePerCore>("Error reading site type");
      }

      SiteReadResult readInSite(isFluid != 0);

      // If solid, there's nothing more to do.
      if (!readInSite.isFluid)
      {
        return readInSite;
      }

      const io::formats::geometry::DisplacementVector& neighbourhood = io::formats::geometry::Get().GetNeighbourhood();
      // Prepare the links array to have enough space.
      readInSite.links.resize(D3Q15::NUMVECTORS - 1);

      // For each link direction...
      for (Direction readDirection = 0; readDirection < neighbourhood.size(); readDirection++)
      {
        // read the type of the intersection and create a link...
        unsigned intersectionType;
        reader.readUnsignedInt(intersectionType);

        LinkReadResult link;
        link.type = (LinkReadResult::IntersectionType) intersectionType;

        // walls have a floating-point distance to the wall...
        if (link.type == LinkReadResult::WALL_INTERSECTION)
        {
          float distance;
          reader.readFloat(distance);
          link.distanceToIntersection = distance;
        }
        // inlets and outlets (which together with none make up the other intersection types)
        // have an iolet id and a distance float...
        else if (link.type != LinkReadResult::NO_INTERSECTION)
        {
          float distance;
          unsigned ioletId;
          reader.readUnsignedInt(ioletId);
          reader.readFloat(distance);

          link.ioletId = ioletId;
          link.distanceToIntersection = distance;
        }

        // Now, attempt to match the direction read from the local neighbourhood to one in the
        // lattice being used for simulation. If a match is found, assign the link to the read
        // site.
        for (Direction usedLatticeDirection = 1; usedLatticeDirection < D3Q15::NUMVECTORS; usedLatticeDirection++)
        {
          if (D3Q15::CX[usedLatticeDirection] == neighbourhood[readDirection].x
              && D3Q15::CY[usedLatticeDirection] == neighbourhood[readDirection].y
              && D3Q15::CZ[usedLatticeDirection] == neighbourhood[readDirection].z)
          {
            // If this link direction is necessary to the lattice in use, keep the link data.
            readInSite.links[usedLatticeDirection - 1] = link;
            break;
          }
        }
      }

      return readInSite;
    }

    proc_t GeometryReader::GetReadingCoreForBlock(site_t blockNumber)
    {
      return proc_t(blockNumber % util::NumericalFunctions::min(READING_GROUP_SIZE, currentCommSize));
    }

    void GeometryReader::ValidateAllReadData()
    {
      if (log::Logger::ShouldDisplay<log::Debug>())
      {
        log::Logger::Log<log::Debug, log::OnePerCore>("Validating the GlobalLatticeData");

        // We check the isFluid property and the link type for each direction
        site_t blockSiteDataLength = readingResult.GetSitesPerBlock() * (1 + D3Q15::NUMVECTORS - 1);

        proc_t* procForSiteRecv = new proc_t[readingResult.GetSitesPerBlock()];
        proc_t * myProcForSite = new proc_t[readingResult.GetSitesPerBlock()];
        unsigned* dummySiteData = new unsigned[blockSiteDataLength];
        unsigned* siteDataRecv = new unsigned[blockSiteDataLength];

        // We also validate that each processor has the same beliefs about each site.
        for (site_t block = 0; block < readingResult.GetBlockCount(); ++block)
        {
          if (readingResult.Blocks[block].Sites.size() == 0)
          {
            for (site_t localSite = 0; localSite < readingResult.GetSitesPerBlock(); ++localSite)
            {
              myProcForSite[localSite] = BIG_NUMBER2;
              dummySiteData[localSite * D3Q15::NUMVECTORS] = std::numeric_limits<unsigned>::max();
              for (Direction direction = 1; direction < D3Q15::NUMVECTORS; ++direction)
              {
                dummySiteData[localSite * D3Q15::NUMVECTORS + direction] = std::numeric_limits<unsigned>::max();
              }
            }
          }
          else
          {
            for (site_t localSite = 0; localSite < readingResult.GetSitesPerBlock(); ++localSite)
            {
              myProcForSite[localSite] = readingResult.Blocks[block].Sites[localSite].targetProcessor;
              dummySiteData[localSite * D3Q15::NUMVECTORS] = readingResult.Blocks[block].Sites[localSite].isFluid;
            
              for (Direction direction = 1; direction < D3Q15::NUMVECTORS; ++direction)
              {
                if (readingResult.Blocks[block].Sites[localSite].isFluid) 
                {
                  dummySiteData[localSite * D3Q15::NUMVECTORS + direction] =
                      readingResult.Blocks[block].Sites[localSite].links[direction - 1].type;
                } 
                else 
                {
                   dummySiteData[localSite * D3Q15::NUMVECTORS + direction] = std::numeric_limits<unsigned>::max();
                }
              }
            }
          }

          // Reduce using a minimum to find the actual processor for each site (ignoring the
          // BIG_NUMBER2 entries).
          MPI_Allreduce(myProcForSite,
                        procForSiteRecv,
                        (int) readingResult.GetSitesPerBlock(),
                        MpiDataType(procForSiteRecv[0]),
                        MPI_MIN,
                        topologyComm);

          MPI_Allreduce(dummySiteData,
                        siteDataRecv,
                        (int) blockSiteDataLength,
                        MpiDataType(dummySiteData[0]),
                        MPI_MIN,
                        topologyComm);

          for (site_t site = 0; site < readingResult.GetSitesPerBlock(); ++site)
          {
            if (procForSiteRecv[site] == ConvertTopologyRankToGlobalRank(topologyRank)
                && (myProcForSite[site] != ConvertTopologyRankToGlobalRank(topologyRank)))
            {
              log::Logger::Log<log::Debug, log::OnePerCore>("Other cores think this core has site %li on block %li but it disagrees.",
                                                            site,
                                                            block);
            }
            else if (myProcForSite[site] != BIG_NUMBER2 && procForSiteRecv[site] != myProcForSite[site])
            {
              log::Logger::Log<log::Debug, log::OnePerCore>("This core thought that core %li has site %li on block %li but others think it's on core %li.",
                                                            myProcForSite[site],
                                                            procForSiteRecv[site],
                                                            site,
                                                            block);
            }

            if (readingResult.Blocks[block].Sites.size() > 0)
            {
              if (dummySiteData[site * D3Q15::NUMVECTORS] != siteDataRecv[site * D3Q15::NUMVECTORS])
              {
                log::Logger::Log<log::Debug, log::OnePerCore>("Different fluid state was found for site %li on block %li. One: %li, Two: %li .",
                                                              site,
                                                              block,
                                                              dummySiteData[site * D3Q15::NUMVECTORS],
                                                              siteDataRecv[site * D3Q15::NUMVECTORS]);
              }

              for (Direction direction = 1; direction < D3Q15::NUMVECTORS; ++direction)
              {
                if (dummySiteData[site * D3Q15::NUMVECTORS + direction]
                    != siteDataRecv[site * D3Q15::NUMVECTORS + direction])
                {
                  log::Logger::Log<log::Debug, log::OnePerCore>("Different link type was found for site %li, link %i on block %li. One: %li, Two: %li .",
                                                                site,
                                                                direction,
                                                                block,
                                                                dummySiteData[site * D3Q15::NUMVECTORS + direction],
                                                                siteDataRecv[site * D3Q15::NUMVECTORS + direction]);
                }
              }

            }
          }
        }

        delete[] procForSiteRecv;
        delete[] myProcForSite;
        delete[] dummySiteData;
        delete[] siteDataRecv;
      }
    }

    /**
     * Compile a list of blocks to be read onto this core, including all the ones we perform
     * LB on, and also any of their neighbouring blocks.
     *
     * NOTE: that the skipping-over of blocks without any fluid sites is dealt with by other
     * code.
     *
     * @param readBlock
     * @param unitForEachBlock
     * @param localRank
     * @param iGlobLatDat
     */
    void GeometryReader::DecideWhichBlocksToRead(bool* readBlock,
                                                 const proc_t* unitForEachBlock,
                                                 const proc_t localRank)
    {
      // Initialise each block to not being read.
      for (site_t block = 0; block < readingResult.GetBlockCount(); ++block)
      {
        readBlock[block] = false;
      }

      // Read a block in if it has fluid sites and is to live on the current processor. Also read
      // in any neighbours with fluid sites.
      for (site_t blockI = 0; blockI < readingResult.blocks.x; ++blockI)
      {
        for (site_t blockJ = 0; blockJ < readingResult.blocks.y; ++blockJ)
        {
          for (site_t blockK = 0; blockK < readingResult.blocks.z; ++blockK)
          {
            site_t lBlockId = readingResult.GetBlockIdFromBlockCoordinates(blockI, blockJ, blockK);

            if (unitForEachBlock[lBlockId] != localRank)
            {
              continue;
            }

            // Read in all neighbouring blocks.
            for (site_t neighI = util::NumericalFunctions::max<site_t>(0, blockI - 1);
                (neighI <= (blockI + 1)) && (neighI < readingResult.blocks.x); ++neighI)
            {
              for (site_t neighJ = util::NumericalFunctions::max<site_t>(0, blockJ - 1);
                  (neighJ <= (blockJ + 1)) && (neighJ < readingResult.blocks.y); ++neighJ)
              {
                for (site_t neighK = util::NumericalFunctions::max<site_t>(0, blockK - 1);
                    (neighK <= (blockK + 1)) && (neighK < readingResult.blocks.z); ++neighK)
                {
                  site_t lNeighId = readingResult.GetBlockIdFromBlockCoordinates(neighI, neighJ, neighK);

                  readBlock[lNeighId] = true;
                }
              }
            }
          }
        }
      }
    }

    /**
     * Get an initial distribution of which block should be on which processor.
     *
     * To make this fast and remove a need to read in the site info about blocks,
     * we assume that neighbouring blocks with any fluid sites on have lattice links
     * between them.
     *
     * NOTE that the old version of this code used to cope with running on multiple machines,
     * by decomposing fluid sites over machines, then decomposing over processors within one machine.
     * To achieve this here, use parmetis's "nparts" parameters to first decompose over machines, then
     * to decompose within machines.
     *
     * @param initialProcForEachBlock Array of length iBlockCount, into which the rank number will be
     * written for each block
     * @param blockCountPerProc Array of length topology size, into which the number of blocks
     * allocated to each processor will be written.
     */
    void GeometryReader::BlockDecomposition()
    {
      site_t* blockCountPerProc = new site_t[topologySize];

      // Count of block sites.
      site_t unvisitedFluidBlockCount = 0;
      for (site_t block = 0; block < readingResult.GetBlockCount(); ++block)
      {
        if (fluidSitesPerBlock[block] != 0)
        {
          ++unvisitedFluidBlockCount;
        }
      }

      // Initialise site count per processor
      for (unsigned int rank = 0; rank < topologySize; ++rank)
      {
        blockCountPerProc[rank] = 0;
      }

      // Divide blocks between the processors.
      DivideBlocks(unvisitedFluidBlockCount,
                   topologySize,
                   blockCountPerProc,
                   procForEachBlock,
                   fluidSitesPerBlock);

      delete[] blockCountPerProc;
    }

    void GeometryReader::ValidateProcForEachBlock()
    {
      if (log::Logger::ShouldDisplay<log::Debug>())
      {
        log::Logger::Log<log::Debug, log::OnePerCore>("Validating procForEachBlock");

        proc_t* procForEachBlockRecv = new proc_t[readingResult.GetBlockCount()];

        MPI_Allreduce(procForEachBlock,
                      procForEachBlockRecv,
                      (int) readingResult.GetBlockCount(),
                      MpiDataType<proc_t>(),
                      MPI_MAX,
                      topologyComm);

        for (site_t block = 0; block < readingResult.GetBlockCount(); ++block)
        {
          if (procForEachBlock[block] != procForEachBlockRecv[block])
          {
            log::Logger::Log<log::Debug, log::OnePerCore>("At least one other proc thought block %li should be on proc %li but we locally had it as %li",
                                                          block,
                                                          procForEachBlock[block],
                                                          procForEachBlockRecv[block]);
          }
        }

        delete[] procForEachBlockRecv;
      }
    }

    /**
     * Get an initial decomposition of the domain macro-blocks over processors (or some other
     * unit of computation, like an entire machine).
     *
     * NOTE: We need the global lattice data and fluid sites per block in order to try to keep
     * contiguous blocks together, and to skip blocks with no fluid sites.
     *
     * @param unassignedBlocks
     * @param blockCount
     * @param unitCount
     * @param blocksOnEachUnit
     * @param unitForEachBlock
     * @param fluidSitesPerBlock
     * @param iGlobLatDat
     */
    void GeometryReader::DivideBlocks(site_t unassignedBlocks,
                                      const proc_t unitCount,
                                      site_t* blocksOnEachUnit,
                                      proc_t* unitForEachBlock,
                                      const site_t* fluidSitesPerBlock)
    {
      // Initialise the unit being assigned to, and the approximate number of blocks
      // required on each unit.
      proc_t currentUnit = 0;

      site_t blocksPerUnit = (site_t) ceil((double) unassignedBlocks / (double) (topologySize));

      // Create an array to monitor whether each block has been assigned yet.
      bool *blockAssigned = new bool[readingResult.GetBlockCount()];

      for (site_t block = 0; block < readingResult.GetBlockCount(); ++block)
      {
        blockAssigned[block] = false;
      }

      std::vector < util::Vector3D<site_t> > currentEdge;
      std::vector < util::Vector3D<site_t> > expandedEdge;

      site_t blockNumber = -1;

      // Domain Decomposition.  Pick a site. Set it to the rank we are
      // looking at. Find its neighbours and put those on the same
      // rank, then find the next-nearest neighbours, etc. until we
      // have a completely joined region, or there are enough fluid
      // sites on the rank.  In the former case, start again at
      // another site. In the latter case, move on to the next rank.
      // Do this until all sites are assigned to a rank. There is a
      // high chance of of all sites on a rank being joined.

      site_t blocksOnCurrentProc = 0;

      // Iterate over all blocks.
      for (site_t blockCoordI = 0; blockCoordI < readingResult.blocks.x; blockCoordI++)
      {
        for (site_t blockCoordJ = 0; blockCoordJ < readingResult.blocks.y; blockCoordJ++)
        {
          for (site_t blockCoordK = 0; blockCoordK < readingResult.blocks.z; blockCoordK++)
          {
            // Block number is the number of the block we're currently on.
            blockNumber++;

            // If the array of proc rank for each site is NULL, we're on an all-solid block.
            // Alternatively, if this block has already been assigned, move on.
            if (fluidSitesPerBlock[blockNumber] == 0)
            {
              unitForEachBlock[blockNumber] = -1;
              continue;
            }
            else if (blockAssigned[blockNumber])
            {
              continue;
            }

            // Assign this block to the current unit.
            blockAssigned[blockNumber] = true;
            unitForEachBlock[blockNumber] = currentUnit;

            ++blocksOnCurrentProc;

            // Record the location of this initial site.
            currentEdge.clear();
            util::Vector3D<site_t> lNew;
            lNew.x = blockCoordI;
            lNew.y = blockCoordJ;
            lNew.z = blockCoordK;
            currentEdge.push_back(lNew);

            // The subdomain can grow.
            bool isRegionGrowing = true;

            // While the region can grow (i.e. it is not bounded by solids or visited
            // sites), and we need more sites on this particular rank.
            while (blocksOnCurrentProc < blocksPerUnit && isRegionGrowing)
            {
              expandedEdge.clear();

              // Sites added to the edge of the mClusters during the iteration.
              isRegionGrowing = Expand(currentEdge,
                                       expandedEdge,
                                       fluidSitesPerBlock,
                                       blockAssigned,
                                       currentUnit,
                                       unitForEachBlock,
                                       blocksOnCurrentProc,
                                       blocksPerUnit);

              // When the new layer of edge sites has been found, swap the buffers for
              // the current and new layers of edge sites.
              currentEdge.swap(expandedEdge);
            }

            // If we have enough sites, we have finished.
            if (blocksOnCurrentProc >= blocksPerUnit)
            {
              blocksOnEachUnit[currentUnit] = blocksOnCurrentProc;

              ++currentUnit;

              unassignedBlocks -= blocksOnCurrentProc;
              blocksPerUnit = (site_t) ceil((double) unassignedBlocks / (double) (unitCount - currentUnit));

              blocksOnCurrentProc = 0;
            }
            // If not, we have to start growing a different region for the same rank:
            // region expansions could get trapped.

          } // Block co-ord k
        } // Block co-ord j
      } // Block co-ord i

      delete[] blockAssigned;
    }

    /**
     * Returns true if the region was expanded.
     *
     * @param edgeBlocks
     * @param iGlobLatDat
     * @param fluidSitesPerBlock
     * @param blockAssigned
     * @param currentUnit
     * @param unitForEachBlock
     * @param blocksOnCurrentProc
     * @param blocksPerUnit
     * @return
     */
    bool GeometryReader::Expand(std::vector<util::Vector3D<site_t> >& edgeBlocks,
                                std::vector<util::Vector3D<site_t> >& expansionBlocks,
                                const site_t* fluidSitesPerBlock,
                                bool* blockAssigned,
                                const proc_t currentUnit,
                                proc_t* unitForEachBlock,
                                site_t &blocksOnCurrentUnit,
                                const site_t blocksPerUnit)
    {
      bool regionExpanded = false;

      // For sites on the edge of the domain (sites_a), deal with the neighbours.
      for (unsigned int index_a = 0; index_a < edgeBlocks.size() && blocksOnCurrentUnit < blocksPerUnit; index_a++)
      {
        util::Vector3D<site_t>& lNew = edgeBlocks[index_a];

        for (unsigned int l = 1; l < D3Q15::NUMVECTORS && blocksOnCurrentUnit < blocksPerUnit; l++)
        {
          // Record neighbour location.
          site_t neighbourI = lNew.x + D3Q15::CX[l];
          site_t neighbourJ = lNew.y + D3Q15::CY[l];
          site_t neighbourK = lNew.z + D3Q15::CZ[l];

          // Move on if neighbour is outside the bounding box.
          if (neighbourI == -1 || neighbourI == readingResult.blocks.x)
            continue;
          if (neighbourJ == -1 || neighbourJ == readingResult.blocks.y)
            continue;
          if (neighbourK == -1 || neighbourK == readingResult.blocks.z)
            continue;

          // Move on if the neighbour is in a block of solids (in which case
          // the pointer to ProcessorRankForEachBlockSite is NULL) or it is solid or has already
          // been assigned to a rank (in which case ProcessorRankForEachBlockSite != -1).  ProcessorRankForEachBlockSite
          // was initialized in lbmReadConfig in io.cc.

          site_t neighBlockId = readingResult.GetBlockIdFromBlockCoordinates(neighbourI, neighbourJ, neighbourK);

          // Don't use this block if it has no fluid sites, or if it has already been assigned to a processor.
          if (fluidSitesPerBlock[neighBlockId] == 0 || blockAssigned[neighBlockId])
          {
            continue;
          }

          // Set the rank for a neighbour and update the fluid site counters.
          blockAssigned[neighBlockId] = true;
          unitForEachBlock[neighBlockId] = currentUnit;
          ++blocksOnCurrentUnit;

          // Neighbour was found, so the region can grow.
          regionExpanded = true;

          // Record the location of the neighbour.
          util::Vector3D<site_t> lNewB(neighbourI, neighbourJ, neighbourK);
          expansionBlocks.push_back(lNewB);
        }
      }

      return regionExpanded;
    }

    void GeometryReader::OptimiseDomainDecomposition(const proc_t* procForEachBlock)
    {
      // Get some arrays that ParMetis needs.
      idx_t* vtxDistribn = new idx_t[topologySize + 1];

      GetSiteDistributionArray(vtxDistribn, procForEachBlock, fluidSitesPerBlock);

      idx_t* firstSiteIndexPerBlock = new idx_t[readingResult.GetBlockCount()];

      GetFirstSiteIndexOnEachBlock(firstSiteIndexPerBlock, vtxDistribn, procForEachBlock, fluidSitesPerBlock);

      idx_t localVertexCount = vtxDistribn[topologyRank + 1] - vtxDistribn[topologyRank];

      idx_t* adjacenciesPerVertex = new idx_t[localVertexCount + 1U];
      idx_t* localAdjacencies;

      GetAdjacencyData(adjacenciesPerVertex,
                       localAdjacencies,
                       localVertexCount,
                       procForEachBlock,
                       firstSiteIndexPerBlock);

      if (log::Logger::ShouldDisplay<log::Debug>())
      {
        ValidateGraphData(vtxDistribn, localVertexCount, adjacenciesPerVertex, localAdjacencies);
      }

      // Call parmetis.
      idx_t* partitionVector = new idx_t[localVertexCount];
      timings[hemelb::reporting::Timers::parmetis].Start();
      CallParmetis(partitionVector, localVertexCount, vtxDistribn, adjacenciesPerVertex, localAdjacencies);
      timings[hemelb::reporting::Timers::parmetis].Stop();
      delete[] localAdjacencies;
      delete[] adjacenciesPerVertex;
      
      // Convert the ParMetis results into a nice format.
      idx_t* allMoves = new idx_t[topologySize];
      
      idx_t* movesList = GetMovesList(allMoves,
                                      firstSiteIndexPerBlock,
                                      procForEachBlock,
                                      fluidSitesPerBlock,
                                      vtxDistribn,
                                      partitionVector);
      
      delete[] firstSiteIndexPerBlock;
      delete[] vtxDistribn;
      delete[] partitionVector;
      timings[hemelb::reporting::Timers::reRead].Start();
      // Reread the blocks based on the ParMetis decomposition.
      RereadBlocks(allMoves, movesList, procForEachBlock);
      timings[hemelb::reporting::Timers::reRead].Stop();
      
      timings[hemelb::reporting::Timers::moves].Start();
      // Implement the decomposition now that we have read the necessary data.
      ImplementMoves(procForEachBlock, allMoves, movesList);
      timings[hemelb::reporting::Timers::moves].Stop();
      delete[] movesList;
      delete[] allMoves;
    }

    // The header section of the config file contains a number of records.
    site_t GeometryReader::GetHeaderLength(site_t blockCount) const
    {
      return io::formats::geometry::HeaderRecordLength * blockCount;
    }

    /**
     * Get the cumulative count of sites on each processor.
     *
     * @param vertexDistribn
     * @param blockCount
     * @param procForEachBlock
     * @param fluidSitesPerBlock
     */
    void GeometryReader::GetSiteDistributionArray(idx_t* vertexDistribn,
                                                  const proc_t* procForEachBlock,
                                                  const site_t* fluidSitesPerBlock) const
    {
      // Firstly, count the sites per processor.
      for (unsigned int rank = 0; rank < (topologySize + 1); ++rank)
      {
        vertexDistribn[rank] = 0;
      }

      for (site_t block = 0; block < readingResult.GetBlockCount(); ++block)
      {
        if (procForEachBlock[block] >= 0)
        {
          vertexDistribn[1 + procForEachBlock[block]] += (idx_t) fluidSitesPerBlock[block];
        }
      }

      // Now make the count cumulative.
      for (unsigned int rank = 0; rank < topologySize; ++rank)
      {
        vertexDistribn[rank + 1] += vertexDistribn[rank];
      }

      // Validate if we're logging.
      if (log::Logger::ShouldDisplay<log::Debug>())
      {
        log::Logger::Log<log::Debug, log::OnePerCore>("Validating the vertex distribution.");

        // vtxDistribn should be the same on all cores.
        idx_t* vtxDistribnRecv = new idx_t[topologySize + 1];

        MPI_Allreduce(vertexDistribn,
                      vtxDistribnRecv,
                      topologySize + 1,
                      MpiDataType(vtxDistribnRecv[0]),
                      MPI_MIN,
                      topologyComm);

        for (unsigned int rank = 0; rank < topologySize + 1; ++rank)
        {
          if (vertexDistribn[rank] != vtxDistribnRecv[rank])
          {
            log::Logger::Log<log::Debug, log::OnePerCore>("vertexDistribn[%i] was %li but at least one other core had it as %li.",
                                                          rank,
                                                          vertexDistribn[rank],
                                                          vtxDistribnRecv[rank]);
          }
        }

        delete[] vtxDistribnRecv;
      }
    }

    void GeometryReader::GetFirstSiteIndexOnEachBlock(idx_t* firstSiteIndexPerBlock,
                                                      const idx_t* vertexDistribution,
                                                      const proc_t* procForEachBlock,
                                                      const site_t* fluidSitesPerBlock) const
    {
      // First calculate the lowest site index on each proc - relatively easy.
      site_t* firstSiteOnProc = new site_t[topologySize];
      for (unsigned int rank = 0; rank < topologySize; ++rank)
      {
        firstSiteOnProc[rank] = vertexDistribution[rank];
      }

      // Now for each block (in ascending order), the smallest site index is the smallest site
      // index on its processor, incremented by the number of sites observed from that processor
      // so far.
      for (site_t block = 0; block < readingResult.GetBlockCount(); ++block)
      {
        proc_t proc = procForEachBlock[block];
        if (proc < 0)
        {
          firstSiteIndexPerBlock[block] = -1;
        }
        else
        {
          firstSiteIndexPerBlock[block] = (idx_t) firstSiteOnProc[proc];
          firstSiteOnProc[proc] += fluidSitesPerBlock[block];
        }
      }

      // Clean up.
      delete[] firstSiteOnProc;

      // Now, if logging debug info, we validate firstSiteIndexPerBlock
      if (log::Logger::ShouldDisplay<log::Debug>())
      {
        log::Logger::Log<log::Debug, log::OnePerCore>("Validating the firstSiteIndexPerBlock values.");

        idx_t* firstSiteIndexPerBlockRecv = new idx_t[readingResult.GetBlockCount()];

        // Reduce finding the maximum across all nodes. Note that we have to use the maximum
        // because some cores will have -1 for a block (indicating that it has no neighbours on
        // that block.
        MPI_Allreduce(firstSiteIndexPerBlock,
                      firstSiteIndexPerBlockRecv,
                      (int) readingResult.GetBlockCount(),
                      MpiDataType(firstSiteIndexPerBlock[0]),
                      MPI_MAX,
                      topologyComm);

        for (site_t block = 0; block < readingResult.GetBlockCount(); ++block)
        {
          if (firstSiteIndexPerBlock[block] >= 0 && firstSiteIndexPerBlock[block] != firstSiteIndexPerBlockRecv[block])
          {
            log::Logger::Log<log::Debug, log::OnePerCore>("This core had the first site index on block %li as %li but at least one other core had it as %li.",
                                                          block,
                                                          firstSiteIndexPerBlock[block],
                                                          firstSiteIndexPerBlockRecv[block]);
          }
        }
      }
    }

    void GeometryReader::GetAdjacencyData(idx_t* adjacenciesPerVertex,
                                          idx_t* &localAdjacencies,
                                          const idx_t localVertexCount,
                                          const proc_t* procForEachBlock,
                                          const idx_t* firstSiteIndexPerBlock) const
    {
      std::vector < idx_t > adj;

      idx_t totalAdjacencies = 0;
      adjacenciesPerVertex[0] = 0;
      idx_t fluidVertex = 0;

      // For each block (counting up by lowest site id)...
      for (site_t blockI = 0; blockI < readingResult.blocks.x; blockI++)
      {
        for (site_t blockJ = 0; blockJ < readingResult.blocks.y; blockJ++)
        {
          for (site_t blockK = 0; blockK < readingResult.blocks.z; blockK++)
          {
            const site_t blockNumber = readingResult.GetBlockIdFromBlockCoordinates(blockI, blockJ, blockK);

            // ... considering only the ones which live on this proc...
            if (procForEachBlock[blockNumber] != topologyRank)
            {
              continue;
            }

            BlockReadResult& blockReadResult = readingResult.Blocks[blockNumber];

            site_t m = -1;

            // ... iterate over sites within the block...
            for (site_t localSiteI = 0; localSiteI < readingResult.blockSize; localSiteI++)
            {
              for (site_t localSiteJ = 0; localSiteJ < readingResult.blockSize; localSiteJ++)
              {
                for (site_t localSiteK = 0; localSiteK < readingResult.blockSize; localSiteK++)
                {
                  ++m;

                  // ... only looking at non-solid sites...
                  if (blockReadResult.Sites[m].targetProcessor == BIG_NUMBER2)
                  {
                    continue;
                  }

                  // ... for each lattice direction...
                  for (unsigned int l = 1; l < D3Q15::NUMVECTORS; l++)
                  {
                    // ... which leads to a valid neighbouring site...
                    site_t neighbourI = blockI * readingResult.blockSize + localSiteI + D3Q15::CX[l];
                    site_t neighbourJ = blockJ * readingResult.blockSize + localSiteJ + D3Q15::CY[l];
                    site_t neighbourK = blockK * readingResult.blockSize + localSiteK + D3Q15::CZ[l];

                    if (neighbourI < 0 || neighbourJ < 0 || neighbourK < 0
                        || neighbourI >= (readingResult.blockSize * readingResult.blocks.x)
                        || neighbourJ >= (readingResult.blockSize * readingResult.blocks.y)
                        || neighbourK >= (readingResult.blockSize * readingResult.blocks.z))
                    {
                      continue;
                    }

                    // ... (that is actually being simulated and not a solid)...
                    site_t neighbourBlockI = neighbourI / readingResult.blockSize;
                    site_t neighbourBlockJ = neighbourJ / readingResult.blockSize;
                    site_t neighbourBlockK = neighbourK / readingResult.blockSize;

                    site_t neighbourSiteI = neighbourI % readingResult.blockSize;
                    site_t neighbourSiteJ = neighbourJ % readingResult.blockSize;
                    site_t neighbourSiteK = neighbourK % readingResult.blockSize;

                    site_t neighbourBlockId = readingResult.GetBlockIdFromBlockCoordinates(neighbourBlockI,
                                                                                           neighbourBlockJ,
                                                                                           neighbourBlockK);

                    BlockReadResult& neighbourBlock = readingResult.Blocks[neighbourBlockId];

                    site_t neighbourSiteId = readingResult.GetSiteIdFromSiteCoordinates(neighbourSiteI,
                                                                                        neighbourSiteJ,
                                                                                        neighbourSiteK);

                    if (neighbourBlock.Sites.size() == 0
                        || neighbourBlock.Sites[neighbourSiteId].targetProcessor == BIG_NUMBER2)
                    {
                      continue;
                    }

                    // Calculate the site's id over the whole geometry,
                    site_t neighGlobalSiteId = firstSiteIndexPerBlock[neighbourBlockId];

                    for (site_t neighSite = 0; neighSite < readingResult.GetSitesPerBlock(); ++neighSite)
                    {
                      if (neighSite == neighbourSiteId)
                      {
                        break;
                      }
                      else if (neighbourBlock.Sites[neighSite].targetProcessor != BIG_NUMBER2)
                      {
                        ++neighGlobalSiteId;
                      }
                    }

                    // then add this to the list of adjacencies.
                    ++totalAdjacencies;
                    adj.push_back((idx_t) neighGlobalSiteId);
                  }

                  // The cumulative count of adjacencies for this vertex is equal to the total
                  // number of adjacencies we've entered.
                  // NOTE: The prefix operator is correct here because
                  // the array has a leading 0 not relating to any site.
                  adjacenciesPerVertex[++fluidVertex] = (idx_t) totalAdjacencies;
                }
              }
            }
          }
        }
      }

      // Perform a debugging test if running at the appropriate log level.
      if (log::Logger::ShouldDisplay<log::Debug>())
      {
        if (fluidVertex != localVertexCount)
        {
          std::cerr << "Encountered a different number of vertices on two different parses: " << fluidVertex << " and "
              << localVertexCount << "\n";
        }
      }

      localAdjacencies = new idx_t[totalAdjacencies];
      for (idx_t adjacency = 0; adjacency < totalAdjacencies; ++adjacency)
      {
        localAdjacencies[adjacency] = adj[adjacency];
      }
    }

    void GeometryReader::CallParmetis(idx_t* partitionVector,
                                      idx_t localVertexCount,
                                      idx_t* vtxDistribn,
                                      idx_t* adjacenciesPerVertex,
                                      idx_t* adjacencies)
    {
      // From the ParMETIS documentation:
      // --------------------------------
      // Processor Pi holds ni consecutive vertices and mi corresponding edges
      //
      // xadj[ni+1] has the cumulative number of adjacencies per vertex (with a leading 0 on each processor)
      // vwgt[ni] has vertex weight coefficients and can be NULL
      // adjncy[mi] has the adjacent vertices for each edge (using a global index, starting at 0)
      // adjwgt[mi] has edge weights and can be NULL
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
      for (idx_t vertexIndex = 0; vertexIndex < localVertexCount; ++vertexIndex)
      {
        partitionVector[vertexIndex] = topologyRank;
      }

      // Weight all vertices evenly.
      idx_t* vertexWeight = new idx_t[localVertexCount];
      for (idx_t vertexIndex = 0; vertexIndex < localVertexCount; ++vertexIndex)
      {
        vertexWeight[vertexIndex] = 1;
      }

      // Set the weights of each partition to be even, and to sum to 1.
      idx_t desiredPartitionSize = topologySize;

      real_t* domainWeights = new real_t[desiredPartitionSize];
      for (idx_t rank = 0; rank < desiredPartitionSize; ++rank)
      {
        domainWeights[rank] = (real_t)(1.0) / ((real_t) desiredPartitionSize);
      }

      // A bunch of values ParMetis needs.
      idx_t noConstraints = 1;
      idx_t weightFlag = 2;
      idx_t numberingFlag = 0;
      idx_t edgesCut = 0;
      idx_t options[4] = { 0, 0, 0, 0 };

      if (log::Logger::ShouldDisplay<log::Debug>())
      {
        // Specify that some options are set and that we should
        // debug everything.
        options[0] = 1;
        options[1] = 1 | 2 | 4 | 8 | 32 | 64;
      }

      real_t tolerance = 1.001F;

      log::Logger::Log<log::Debug, log::OnePerCore>("Calling ParMetis");
      ParMETIS_V3_PartKway(vtxDistribn,
                           adjacenciesPerVertex,
                           adjacencies,
                           vertexWeight,
                           NULL,
                           &weightFlag,
                           &numberingFlag,
                           &noConstraints,
                           &desiredPartitionSize,
                           domainWeights,
                           &tolerance,
                           options,
                           &edgesCut,
                           partitionVector,
                           &topologyComm);
      log::Logger::Log<log::Debug, log::OnePerCore>("ParMetis returned.");

      delete[] domainWeights;
      delete[] vertexWeight;
    }

    void GeometryReader::ValidateGraphData(idx_t* vtxDistribn,
                                           idx_t localVertexCount,
                                           idx_t* adjacenciesPerVertex,
                                           idx_t* adjacencies)
    {
      // If we're using debugging logs, check that the arguments are consistent across all cores.
      // To verify: vtxDistribn, adjacenciesPerVertex, adjacencies
      if (log::Logger::ShouldDisplay<log::Debug>())
      {
        log::Logger::Log<log::Debug, log::OnePerCore>("Validating the graph adjacency structure");

        // Create an array of lists to store all of this node's adjacencies, arranged by the
        // proc the adjacent vertex is on.
        std::multimap < idx_t, idx_t > *adjByNeighProc = new std::multimap<idx_t, idx_t>[topologySize];
        for (proc_t proc = 0; proc < (proc_t) topologySize; ++proc)
        {
          adjByNeighProc[proc] = std::multimap<idx_t, idx_t>();
        }

        // The adjacency data should correspond across all cores.
        for (idx_t index = 0; index < localVertexCount; ++index)
        {
          idx_t vertex = vtxDistribn[topologyRank] + index;

          // Iterate over each adjacency (of each vertex).
          for (idx_t adjNumber = 0; adjNumber < (adjacenciesPerVertex[index + 1] - adjacenciesPerVertex[index]);
              ++adjNumber)
          {
            idx_t adjacentVertex = adjacencies[adjacenciesPerVertex[index] + adjNumber];
            proc_t adjacentProc = -1;

            // Calculate the proc of the neighbouring vertex.
            for (proc_t proc = 0; proc < (proc_t) topologySize; ++proc)
            {
              if (vtxDistribn[proc] <= adjacentVertex && vtxDistribn[proc + 1] > adjacentVertex)
              {
                adjacentProc = proc;
                break;
              }
            }

            // If it doesn't appear to belong to any proc, something's wrong.
            if (adjacentProc == -1)
            {
              log::Logger::Log<log::Debug, log::OnePerCore>("The vertex %li has a neighbour %li which doesn\'t appear to live on any processor.",
                                                            vertex,
                                                            adjacentVertex);
              continue;
            }

            // Store the data if it does belong to a proc.
            adjByNeighProc[adjacentProc].insert(std::pair<idx_t, idx_t>(adjacentVertex, vertex));
          }
        }

        // Create variables for the neighbour data to go into.
        idx_t* counts = new idx_t[topologySize];
        idx_t** data = new idx_t*[topologySize];
        MPI_Request* requests = new MPI_Request[2 * topologySize];

        log::Logger::Log<log::Debug, log::OnePerCore>("Validating neighbour data");

        // Now spread and compare the adjacency information. Larger ranks send data to smaller
        // ranks which receive the data and compare it.
        for (proc_t neigh = 0; neigh < (proc_t) topologySize; ++neigh)
        {
          data[neigh] = NULL;
          requests[2 * neigh] = MPI_REQUEST_NULL;
          requests[2 * neigh + 1] = MPI_REQUEST_NULL;

          if (neigh < topologyRank)
          {
            // Send the array length.
            counts[neigh] = 2 * adjByNeighProc[neigh].size();
            MPI_Isend(&counts[neigh], 1, MpiDataType(counts[0]), neigh, 42, topologyComm, &requests[2 * neigh]);

            // Create a sendable array (std::lists aren't organised in a sendable format).
            data[neigh] = new idx_t[counts[neigh]];

            unsigned int adjacencyIndex = 0;

            for (std::multimap<idx_t, idx_t>::iterator it = adjByNeighProc[neigh].begin();
                it != adjByNeighProc[neigh].end(); ++it)

            {
              data[neigh][2 * adjacencyIndex] = it->first;
              data[neigh][2 * adjacencyIndex + 1] = it->second;
              ++adjacencyIndex;
            }

            // Send the data to the neighbour.
            MPI_Isend(data[neigh],
                      (int) counts[neigh],
                      MpiDataType<idx_t>(),
                      neigh,
                      43,
                      topologyComm,
                      &requests[2 * neigh + 1]);

            // Sending arrays don't perform comparison.
            continue;
          }

          // If this is a greater rank number than the neighbour, receive the data.
          if (neigh > topologyRank)
          {
            MPI_Irecv(&counts[neigh], 1, MpiDataType(counts[neigh]), neigh, 42, topologyComm, &requests[2 * neigh]);

            MPI_Wait(&requests[2 * neigh], MPI_STATUS_IGNORE);

            data[neigh] = new idx_t[counts[neigh]];

            MPI_Irecv(data[neigh],
                      (int) counts[neigh],
                      MpiDataType<idx_t>(),
                      neigh,
                      43,
                      topologyComm,
                      &requests[2 * neigh + 1]);

            MPI_Wait(&requests[2 * neigh] + 1, MPI_STATUS_IGNORE);
          }
          // Neigh == mTopologyRank, i.e. neighbouring vertices on the same proc
          // Duplicate the data.

          else
          {
            counts[neigh] = 2 * adjByNeighProc[neigh].size();
            data[neigh] = new idx_t[counts[neigh]];

            int adjacencyIndex = 0;
            for (std::multimap<idx_t, idx_t>::iterator it = adjByNeighProc[neigh].begin();
                it != adjByNeighProc[neigh].end(); ++it)

            {
              data[neigh][2 * adjacencyIndex] = it->first;
              data[neigh][2 * adjacencyIndex + 1] = it->second;
              ++adjacencyIndex;
            }
          }

          // Now we compare. First go through the received data which is ordered as (adjacent
          // vertex, vertex) wrt the neighbouring proc.
          for (idx_t ii = 0; ii < counts[neigh]; ii += 2)
          {
            bool found = false;

            // Go through each neighbour we know about on this proc, and check whether it
            // matches the current received neighbour-data.
            for (std::multimap<idx_t, idx_t>::iterator it = adjByNeighProc[neigh].find(data[neigh][ii + 1]);
                it != adjByNeighProc[neigh].end(); ++it)
            {
              idx_t recvAdj = it->first;
              idx_t recvAdj2 = it->second;

              if (data[neigh][ii] == recvAdj2 && data[neigh][ii + 1] == recvAdj)
              {
                adjByNeighProc[neigh].erase(it);
                found = true;
                break;
              }
            }

            // No neighbour data on this proc matched the data received.
            if (!found)
            {
              log::Logger::Log<log::Debug, log::OnePerCore>("Neighbour proc %i had adjacency (%li,%li) that wasn't present on this processor.",
                                                            neigh,
                                                            data[neigh][ii],
                                                            data[neigh][ii + 1]);
            }
          }

          // The local store of adjacencies should now be empty, if there was complete matching.
          std::multimap<idx_t, idx_t>::iterator it = adjByNeighProc[neigh].begin();
          while (it != adjByNeighProc[neigh].end())
          {
            idx_t adj1 = it->first;
            idx_t adj2 = it->second;
            ++it;

            // We had neighbour-data on this proc that didn't match that received.
            log::Logger::Log<log::Debug, log::OnePerCore>("The local processor has adjacency (%li,%li) that isn't present on neighbouring processor %i.",
                                                          adj1,
                                                          adj2,
                                                          neigh);

          }
        }

        // Wait for everything to complete before deallocating.
        MPI_Waitall((int) topologySize, requests, MPI_STATUSES_IGNORE);

        delete[] counts;
        delete[] requests;

        for (proc_t prc = 0; prc < (proc_t) topologySize; ++prc)
        {
          if (data[prc] != NULL)
          {
            delete[] data[prc];
          }
        }
        delete[] data;
        delete[] adjByNeighProc;
      }
    }

    /**
     * Returns a list of the fluid sites to be moved.
     *
     * NOTE: This function's return value is a dynamically-allocated array of all the moves to be
     * performed, ordered by (origin processor [with a count described by the content of the first
     * parameter], site id on the origin processor). The contents of the array are contiguous
     * triplets of ints: (block id, site id on block, destination rank).
     *
     * @param movesFromEachProc
     * @param vtxDistribn
     * @param partitionVector
     * @return
     */
    idx_t* GeometryReader::GetMovesList(idx_t* movesFromEachProc,
                                        const idx_t* firstSiteIndexPerBlock,
                                        const proc_t* procForEachBlock,
                                        const site_t* fluidSitesPerBlock,
                                        const idx_t* vtxDistribn,
                                        const idx_t* partitionVector)
    {
      // Right. Let's count how many sites we're going to have to move. Count the local number of
      // sites to be moved, and collect the site id and the destination processor.
      std::vector < idx_t > moveData;

      const idx_t myLowest = vtxDistribn[topologyRank];
      const idx_t myHighest = vtxDistribn[topologyRank + 1] - 1;

      // For each local fluid site...
      for (idx_t ii = 0; ii <= (myHighest - myLowest); ++ii)
      {
        // ... if it's going elsewhere...
        if (partitionVector[ii] != topologyRank)
        {
          // ... get it's id on the local processor...
          idx_t localFluidSiteId = myLowest + ii;

          // ... find out which block it's on, by going through all blocks until we find one
          // with firstIdOnBlock <= localFluidSiteId < (firstIdOnBlock + sitesOnBlock)...
          idx_t fluidSiteBlock = 0;

          while ( (procForEachBlock[fluidSiteBlock] < 0) || (firstSiteIndexPerBlock[fluidSiteBlock] > localFluidSiteId)
              || ( (firstSiteIndexPerBlock[fluidSiteBlock] + fluidSitesPerBlock[fluidSiteBlock]) <= localFluidSiteId))
          {
            fluidSiteBlock++;
          }

          // ... and find its site id within that block. Start by working out how many fluid sites
          // we have to pass before we arrive at the fluid site we're after...
          idx_t fluidSitesToPass = localFluidSiteId - firstSiteIndexPerBlock[fluidSiteBlock];
          idx_t siteIndex = 0;

          while (true)
          {
            // ... then keep going through the sites on the block until we've passed as many fluid
            // sites as we need to.
            if (readingResult.Blocks[fluidSiteBlock].Sites[siteIndex].targetProcessor != BIG_NUMBER2)
            {
              fluidSitesToPass--;
            }
            if (fluidSitesToPass < 0)
            {
              break;
            }
            siteIndex++;
          }

          // The above code could go wrong, so in debug logging mode, we do some extra tests.
          if (log::Logger::ShouldDisplay<log::Debug>())
          {
            // If we've ended up on an impossible block, or one that doesn't live on this rank,
            // inform the user.
            if (fluidSiteBlock >= readingResult.GetBlockCount() || procForEachBlock[fluidSiteBlock] != topologyRank)
            {
              log::Logger::Log<log::Debug, log::OnePerCore>("Partition element %i wrongly assigned to block %u of %i (block on processor %i)",
                                                            ii,
                                                            fluidSiteBlock,
                                                            readingResult.GetBlockCount(),
                                                            procForEachBlock[fluidSiteBlock]);
            }

            // Similarly, if we've ended up with an impossible site index, or a solid site,
            // print an error message.
            if (siteIndex >= readingResult.GetSitesPerBlock()
                || readingResult.Blocks[fluidSiteBlock].Sites[siteIndex].targetProcessor == BIG_NUMBER2)
            {
              log::Logger::Log<log::Debug, log::OnePerCore>("Partition element %i wrongly assigned to site %u of %i (block %i%s)",
                                                            ii,
                                                            siteIndex,
                                                            fluidSitesPerBlock[fluidSiteBlock],
                                                            fluidSiteBlock,
                                                            readingResult.Blocks[fluidSiteBlock].Sites[siteIndex].targetProcessor
                                                                == BIG_NUMBER2 ?
                                                              " and site is solid" :
                                                              "");
            }
          }

          // Add the block, site and destination rank to our move list.
          moveData.push_back(fluidSiteBlock);
          moveData.push_back(siteIndex);
          moveData.push_back(partitionVector[ii]);
        }
      }

      // Spread the move-count data around, so all processes now how many moves each process is
      // doing.
      idx_t moves = moveData.size() / 3;
      MPI_Allgather(&moves,
                    1,
                    MpiDataType(moves),
                    movesFromEachProc,
                    1,
                    MpiDataType(movesFromEachProc[0]),
                    topologyComm);

      // Count the total moves.
      idx_t totalMoves = 0;

      for (unsigned int ii = 0; ii < topologySize; ++ii)
      {
        totalMoves += movesFromEachProc[ii];
      }

      // Now share all the lists of moves - create a MPI type...
      MPI_Datatype moveType;
      MPI_Type_contiguous(3, MpiDataType<idx_t>(), &moveType);
      MPI_Type_commit(&moveType);

      // ... create a destination array...
      idx_t* movesList = new idx_t[3 * totalMoves];

      // ... create an array of offsets into the destination array for each rank...
      int * offsets = new int[topologySize];

      offsets[0] = 0;
      for (unsigned int ii = 1; ii < topologySize; ++ii)
      {
        offsets[ii] = (int) (offsets[ii - 1] + movesFromEachProc[ii - 1]);
      }

      {
        int* procMovesInt = new int[topologySize];

        for (proc_t ii = 0; ii < (proc_t) topologySize; ++ii)
        {
          procMovesInt[ii] = (int) movesFromEachProc[ii];
        }

        // ... use MPI to gather the data...
        MPI_Allgatherv(&moveData[0], (int) moves, moveType, movesList, procMovesInt, offsets, moveType, topologyComm);

        delete[] procMovesInt;

      }

      // ... clean up...
      MPI_Type_free(&moveType);
      delete[] offsets;

      // ... and return the list of moves.
      return movesList;
    }

    void GeometryReader::RereadBlocks(const idx_t* movesPerProc,
                                      const idx_t* movesList,
                                      const int* procForEachBlock)
    {
      // Initialise the array (of which proc each block belongs to) to what it was before.
      int* newProcForEachBlock = new int[readingResult.GetBlockCount()];

      for (site_t blockNumber = 0; blockNumber < readingResult.GetBlockCount(); ++blockNumber)
      {
        newProcForEachBlock[blockNumber] = procForEachBlock[blockNumber];
      }

      // Set the proc for each block to be the current proc whenever a site on that block is
      // going to be moved to the current proc.
      idx_t moveIndex = 0;

      for (unsigned int fromProc = 0; fromProc < topologySize; ++fromProc)
      {
        for (idx_t moveNumber = 0; moveNumber < movesPerProc[fromProc]; ++moveNumber)
        {
          idx_t block = movesList[3 * moveIndex];
          idx_t toProc = movesList[3 * moveIndex + 2];
          ++moveIndex;

          if (toProc == (idx_t) topologyRank)
          {
            newProcForEachBlock[block] = topologyRank;
          }
        }
      }

      // Reread the blocks into the GlobalLatticeData now.
      ReadInLocalBlocks(newProcForEachBlock, topologyRank);

      // Clean up.
      delete[] newProcForEachBlock;
    }

    void GeometryReader::ImplementMoves(const proc_t* procForEachBlock,
                                        const idx_t* movesFromEachProc,
                                        const idx_t* movesList) const
    {
      // First all, set the proc rank for each site to what it originally was before
      // domain decomposition optimisation. Go through each block...
      for (site_t block = 0; block < readingResult.GetBlockCount(); ++block)
      {
        // If this proc has owned a fluid site on this block either before or after optimisation,
        // the following will be non-null.
        if (readingResult.Blocks[block].Sites.size() > 0)
        {
          // Get the original proc for that block.
          proc_t originalProc = procForEachBlock[block];

          // For each site on that block...
          for (site_t siteIndex = 0; siteIndex < readingResult.GetSitesPerBlock(); ++siteIndex)
          {
            // ... if the site is non-solid...
            if (readingResult.Blocks[block].Sites[siteIndex].targetProcessor != BIG_NUMBER2)
            {
              // ... set its rank to be the rank it had before optimisation.
              readingResult.Blocks[block].Sites[siteIndex].targetProcessor =
                  ConvertTopologyRankToGlobalRank(originalProc);
            }
          }
        }
      }

      // Now implement the moves suggested by parmetis.
      idx_t moveIndex = 0;

      // For each source proc, go through as many moves as it had.
      for (unsigned int fromProc = 0; fromProc < topologySize; ++fromProc)
      {
        for (idx_t moveNumber = 0; moveNumber < movesFromEachProc[fromProc]; ++moveNumber)
        {
          // For each move, get the block, site and destination proc.
          idx_t block = movesList[3 * moveIndex];
          idx_t site = movesList[3 * moveIndex + 1];
          idx_t toProc = movesList[3 * moveIndex + 2];

          // Only implement the move if we have read that block's data.
          if (readingResult.Blocks[block].Sites.size() > 0)
          {
            // Some logging code - the unmodified rank for each move's site should equal
            // lFromProc.
            if (log::Logger::ShouldDisplay<log::Debug>())
            {
              if (readingResult.Blocks[block].Sites[site].targetProcessor
                  != ConvertTopologyRankToGlobalRank((proc_t) fromProc))
              {
                log::Logger::Log<log::Debug, log::OnePerCore>("Block %ld, site %ld from move %u was originally on proc %i, not proc %u.",
                                                              block,
                                                              site,
                                                              moveIndex,
                                                              readingResult.Blocks[block].Sites[site].targetProcessor,
                                                              fromProc);
              }
            }

            // Implement the move.
            readingResult.Blocks[block].Sites[site].targetProcessor = ConvertTopologyRankToGlobalRank((proc_t) toProc);
          }

          ++moveIndex;
        }
      }
    }

    proc_t GeometryReader::ConvertTopologyRankToGlobalRank(proc_t topologyRankIn) const
    {
      // If the global rank is not equal to the topology rank, we are not using rank 0 for
      // LBM.
      return (topology::NetworkTopology::Instance()->GetLocalRank() == topologyRank) ?
        topologyRankIn :
        (topologyRankIn + 1);
    }

  }
}
