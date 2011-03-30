#include <math.h>
#include <parmetis.h>

#include "TopologyReader.h"
#include "io/XdrMemReader.h"
#include "lb/GlobalLatticeData.h"
#include "util/utilityFunctions.h"

namespace hemelb
{
  namespace topology
  {

    // TODO The memory usage here is not yet optimal. We shouldn't allocate any memory for blocks
    // that aren't on this processor or a neighbour.

    // TODO This file is generally ugly. The first step to making it nicer would probably be to
    // integrate it with the LocalLatDat and GlobalLatDat objects (because we really only want
    // LocalLatticeData), and the functions in Net which initialise them. Once the interface
    // to this object is nice and clean, we can tidy up the code here.

    TopologyReader::TopologyReader()
    {
      MPI_Comm_group(MPI_COMM_WORLD, &mTopologyGroup);

      MPI_Group lNewGroup;
      int lExclusions[1] = { 0 };
      MPI_Group_excl(mTopologyGroup, 1, lExclusions, &lNewGroup);
      mTopologyGroup = lNewGroup;

      MPI_Comm_create(MPI_COMM_WORLD, mTopologyGroup, &mTopologyComm);

      MPI_Comm_rank(MPI_COMM_WORLD, &mGlobalRank);

      if (mGlobalRank != 0)
      {
        MPI_Comm_rank(mTopologyComm, &mTopologyRank);
        MPI_Comm_size(mTopologyComm, &mTopologySize);
      }

      // Reduce the size
      --mTopologySize;
    }

    TopologyReader::~TopologyReader()
    {
      MPI_Group_free(&mTopologyGroup);
      if (mGlobalRank != 0)
      {
        MPI_Comm_free(&mTopologyComm);
      }
    }

    /**
     * Read in the section at the beginning of the config file.
     *
     * // PREAMBLE
     * uint stress_type = load_uint();
     * uint blocks[3]; // blocks in x,y,z
     * for (int i=0; i<3; ++i)
     *        blocks[i] = load_uint();
     * // number of sites along 1 edge of a block
     * uint block_size = load_uint();
     *
     * uint nBlocks = blocks[0] * blocks[1] * blocks[2];
     */
    void TopologyReader::ReadPreamble(MPI_File xiFile,
                                      hemelb::lb::LbmParameters * bParams,
                                      hemelb::lb::GlobalLatticeData* bGlobalLatticeData)
    {
      std::string lMode = "native";

      MPI_File_set_view(xiFile, 0, MPI_BYTE, MPI_BYTE, &lMode[0], MPI_INFO_NULL);

      char lPreambleBuffer[PreambleBytes];

      MPI_Status lStatus;

      MPI_File_read_all(xiFile, lPreambleBuffer, PreambleBytes, MPI_BYTE, &lStatus);

      hemelb::io::XdrReader preambleReader = hemelb::io::XdrMemReader(lPreambleBuffer,
                                                                      PreambleBytes);

      // Variables we'll read.
      unsigned int stressType, blocksX, blocksY, blocksZ, blockSize;
      double voxelSize, siteZeroWorldPosition[3];

      preambleReader.readUnsignedInt(stressType);
      preambleReader.readUnsignedInt(blocksX);
      preambleReader.readUnsignedInt(blocksY);
      preambleReader.readUnsignedInt(blocksZ);
      preambleReader.readUnsignedInt(blockSize);
      preambleReader.readDouble(voxelSize);
      for (unsigned int i = 0; i < 3; ++i)
        preambleReader.readDouble(siteZeroWorldPosition[0]);

      bParams->StressType = (lb::StressTypes) stressType;

      bGlobalLatticeData->SetBasicDetails(blocksX, blocksY, blocksZ, blockSize);
    }

    /**
     * Read the header section, with minimal information about each block.
     *
     * // HEADER
     * uint nSites[nBlocks];
     * uint nBytes[nBlocks];
     *
     * for (int i = 0; i < nBlocks; ++i) {
     *        nSites[i] = load_uint(); // sites in block
     *        nBytes[i] = load_uint(); // length, in bytes, of the block's record in this file
     * }
     */
    void TopologyReader::ReadHeader(MPI_File xiFile,
                                    unsigned int iBlockCount,
                                    unsigned int* sitesInEachBlock,
                                    unsigned int* bytesUsedByBlockInDataFile)
    {
      // The header section of the config file contains two adjacent unsigned ints for each block.
      // The first is the number of sites in that block, the second is the number of bytes used by
      // the block in the data file.
      unsigned int headerByteCount = 2 * 4 * iBlockCount;

      char* lHeaderBuffer = new char[headerByteCount];

      MPI_Status lStatus;

      MPI_File_read_all(xiFile, lHeaderBuffer, headerByteCount, MPI_BYTE, &lStatus);

      hemelb::io::XdrReader preambleReader = hemelb::io::XdrMemReader(lHeaderBuffer,
                                                                      headerByteCount);

      for (unsigned int ii = 0; ii < iBlockCount; ii++)
      {
        preambleReader.readUnsignedInt(sitesInEachBlock[ii]);
        preambleReader.readUnsignedInt(bytesUsedByBlockInDataFile[ii]);
      }

      delete[] lHeaderBuffer;
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
     * @param iBlockCount
     * @param sitesInEachBlock
     * @param initialProcForEachBlock
     */
    void TopologyReader::BlockDecomposition(const unsigned int iBlockCount,
                                            const unsigned int iProcCount,
                                            const bool reservedSteeringCore,
                                            const hemelb::lb::GlobalLatticeData* iGlobLatDat,
                                            const unsigned int* fluidSitePerBlock,
                                            int* initialProcForEachBlock,
                                            unsigned int* blockCountPerProc)
    {
      // Count of block sites.
      unsigned int lUnvisitedFluidBlockCount = 0;
      for (unsigned int ii = 0; ii < iBlockCount; ++ii)
      {
        if (fluidSitePerBlock[ii] != 0)
        {
          ++lUnvisitedFluidBlockCount;
        }
      }

      // Initialise site count per processor
      for (unsigned int ii = 0; ii < iProcCount; ++ii)
      {
        blockCountPerProc[ii] = 0;
      }

      // Fluid sites per rank.
      unsigned int blocksPerUnit = (unsigned int) ceil((double) lUnvisitedFluidBlockCount
          / (double) iProcCount);

      //Rank we're looking at.
      unsigned int firstProc = 0;

      // If we're steering with more than one processor, save one processor for doing that.
      if (reservedSteeringCore && iProcCount != 1)
      {
        blocksPerUnit = (unsigned int) ceil((double) lUnvisitedFluidBlockCount
            / (double) (iProcCount - 1));
        firstProc = 1;
        blockCountPerProc[0] = 0;
      }

      // Divide blocks between the processors.
      DivideBlocks(firstProc, blocksPerUnit, lUnvisitedFluidBlockCount, iBlockCount, iProcCount,
                   blockCountPerProc, initialProcForEachBlock, fluidSitePerBlock, iGlobLatDat);
    }

    void TopologyReader::LoadAndDecompose(lb::GlobalLatticeData* bGlobLatDat,
                                          int &totalFluidSites,
                                          unsigned int siteMins[3],
                                          unsigned int siteMaxes[3],
                                          bool iReserveSteeringCore,
                                          NetworkTopology* bNetTop,
                                          lb::LbmParameters* bLbmParams,
                                          SimConfig* bSimConfig,
                                          double* oReadTime,
                                          double* oDecomposeTime)
    {
      double lStart = MPI_Wtime();

      MPI_File lFile;
      MPI_Info lFileInfo;
      int lError;

      // Open the file using the MPI parallel I/O interface at the path
      // given, in read-only mode.
      MPI_Info_create(&lFileInfo);

      // Create hints about how we'll read the file. See Chapter 13, page 400 of the MPI 2.2 spec.
      //
      MPI_Info_set(lFileInfo, "access_style", "read_once,sequential");
      MPI_Info_set(lFileInfo, "collective_buffering", "true");

      lError = MPI_File_open(MPI_COMM_WORLD, &bSimConfig->DataFilePath[0], MPI_MODE_RDONLY,
                             lFileInfo, &lFile);

      if (lError != 0)
      {
        fprintf(stderr, "Unable to open file %s [rank %i], exiting\n",
                bSimConfig->DataFilePath.c_str(), mGlobalRank);
        fflush(0x0);
        exit(0x0);
      }
      else
      {
        fprintf(stderr, "Opened config file %s [rank %i]\n", bSimConfig->DataFilePath.c_str(),
                mGlobalRank);
      }
      fflush(NULL);

      ReadPreamble(lFile, bLbmParams, bGlobLatDat);

      unsigned int* sitesPerBlock = new unsigned int[bGlobLatDat->GetBlockCount()];
      unsigned int* bytesPerBlock = new unsigned int[bGlobLatDat->GetBlockCount()];

      ReadHeader(lFile, bGlobLatDat->GetBlockCount(), sitesPerBlock, bytesPerBlock);

      int* procForEachBlock = new int[bGlobLatDat->GetBlockCount()];
      unsigned int* blockCountForEachProc = new unsigned int[bNetTop->GetProcessorCount()];
      BlockDecomposition(bGlobLatDat->GetBlockCount(), bNetTop->GetProcessorCount(),
                         iReserveSteeringCore, bGlobLatDat, sitesPerBlock, procForEachBlock,
                         blockCountForEachProc);

      ReadInLocalBlocks(lFile, bytesPerBlock, procForEachBlock, bNetTop->GetLocalRank(),
                        bGlobLatDat);

      MPI_File_close(&lFile);
      MPI_Info_free(&lFileInfo);

      double lMiddle = MPI_Wtime();

      int preOptimisationSites = 0;
      for (unsigned int ii = 0; ii < bGlobLatDat->GetBlockCount(); ii++)
      {
        if (procForEachBlock[ii] == (int) bNetTop->GetLocalRank())
        {
          preOptimisationSites += sitesPerBlock[ii];
        }
      }

      OptimiseDomainDecomposition(sitesPerBlock, procForEachBlock, blockCountForEachProc, bNetTop,
                                  bSimConfig, bLbmParams, bGlobLatDat);

      bNetTop->FluidSitesOnEachProcessor = new int[bNetTop->GetProcessorCount()];

      int localFluidSites = 0;

      for (unsigned int lBlock = 0; lBlock < bGlobLatDat->GetBlockCount(); ++lBlock)
      {
        if (bGlobLatDat->Blocks[lBlock].ProcessorRankForEachBlockSite != NULL)
        {
          for (unsigned int lSiteIndex = 0; lSiteIndex < bGlobLatDat->SitesPerBlockVolumeUnit; ++lSiteIndex)
          {
            if (bGlobLatDat->Blocks[lBlock].ProcessorRankForEachBlockSite[lSiteIndex]
                == (int) bNetTop->GetLocalRank())
            {
              ++localFluidSites;
            }
          }
        }
      }

      MPI_Allgather(&localFluidSites, 1, MPI_INT, bNetTop->FluidSitesOnEachProcessor, 1, MPI_INT,
                    MPI_COMM_WORLD);

      //TODO this is a total hack just for now.
      unsigned int localMins[3];
      unsigned int localMaxes[3];
      localMins[0] = UINT_MAX;
      localMins[1] = UINT_MAX;
      localMins[2] = UINT_MAX;
      localMaxes[0] = 0;
      localMaxes[1] = 0;
      localMaxes[2] = 0;

      for (unsigned int siteI = 0; siteI < bGlobLatDat->GetXSiteCount(); ++siteI)
      {
        for (unsigned int siteJ = 0; siteJ < bGlobLatDat->GetYSiteCount(); ++siteJ)
        {
          for (unsigned int siteK = 0; siteK < bGlobLatDat->GetZSiteCount(); ++siteK)
          {
            int* procId = bGlobLatDat->GetProcIdFromGlobalCoords(siteI, siteJ, siteK);
            if (procId == NULL || *procId != (int) bNetTop->GetLocalRank())
            {
              continue;
            }
            else
            {
              localMins[0] = hemelb::util::NumericalFunctions::min<unsigned int>(localMins[0],
                                                                                 siteI);
              localMins[1] = hemelb::util::NumericalFunctions::min<unsigned int>(localMins[1],
                                                                                 siteJ);
              localMins[2] = hemelb::util::NumericalFunctions::min<unsigned int>(localMins[2],
                                                                                 siteK);
              localMaxes[0] = hemelb::util::NumericalFunctions::max<unsigned int>(localMaxes[0],
                                                                                  siteI);
              localMaxes[1] = hemelb::util::NumericalFunctions::max<unsigned int>(localMaxes[1],
                                                                                  siteJ);
              localMaxes[2] = hemelb::util::NumericalFunctions::max<unsigned int>(localMaxes[2],
                                                                                  siteK);
            }
          }
        }
      }

      MPI_Allreduce(localMins, siteMins, 3, MPI_UNSIGNED, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(localMaxes, siteMaxes, 3, MPI_UNSIGNED, MPI_MAX, MPI_COMM_WORLD);

      //TODO this is a total hack just for now.
      totalFluidSites = 0;
      for (unsigned int ii = 0; ii < bGlobLatDat->GetBlockCount(); ++ii)
      {
        totalFluidSites += sitesPerBlock[ii];
      }

      double lEnd = MPI_Wtime();

      *oReadTime = lMiddle - lStart;
      *oDecomposeTime = lEnd - lMiddle;
    }

    /**
     * Read in the necessary blocks from the file.
     *
     * NOTE: we used to set the values of the site_mins and site_maxes (min and max site
     * ids in each of the x, y and z directions) in this function. If this is still
     * needed it probably has to go in the header of the file.
     *
     * BODY:
     * Block *block = new Block[nBlocks];
     * sitesPerBlock = block_size*block_size*block_size;
     *
     * for (int i = 0; i < blocks[0] * blocks[1] * blocks[2]; ++i) {
     *   if (nSites[i] == 0)
     *     // nothing
     *     continue;
     *
     *   block[i]->sites = new Site[sitesPerBlock];
     *   for (int j = 0; j < sitesPerBlock; ++j) {
     *     block[i]->sites[j]->config = load_uint();
     *     if (block[i]->sites[j]->IsAdjacentInletOrOutlet()) {
     *       for (k=0; k<3; ++k)
     *         block[i]->sites[j]->boundaryNormal[k] = load_double();
     *       block[i]->sites[j]->boundaryDistance = load_double();
     *     }
     *
     *     if (block[i]->sites[j]->IsAdjacentToWall()) {
     *       for (k=0; k<3; ++k)
     *         block[i]->sites[j]->wallNormal[k] = load_double();
     *       block[i]->sites[j]->wallDistance = load_double();
     *     }
     *
     *     for (k=0; k<14; ++k)
     *       block[i]->sites[j]->cutDistance[k] = load_double();
     *   }
     * }
     */
    void TopologyReader::ReadInLocalBlocks(MPI_File iFile,
                                           const unsigned int* bytesPerBlock,
                                           const int* unitForEachBlock,
                                           const unsigned int localRank,
                                           const lb::GlobalLatticeData* iGlobLatDat)
    {
      for (unsigned int ii = 0; ii < iGlobLatDat->GetBlockCount(); ++ii)
      {
        if (iGlobLatDat->Blocks[ii].ProcessorRankForEachBlockSite != NULL)
        {
          delete[] iGlobLatDat->Blocks[ii].ProcessorRankForEachBlockSite;
          iGlobLatDat->Blocks[ii].ProcessorRankForEachBlockSite = NULL;
        }
        if (iGlobLatDat->Blocks[ii].site_data != NULL)
        {
          delete[] iGlobLatDat->Blocks[ii].site_data;
          iGlobLatDat->Blocks[ii].site_data = NULL;
        }
        if (iGlobLatDat->Blocks[ii].wall_data != NULL)
        {
          delete[] iGlobLatDat->Blocks[ii].wall_data;
          iGlobLatDat->Blocks[ii].wall_data = NULL;
        }
      }

      // Create a list of which blocks to read in.
      bool* readBlock = new bool[iGlobLatDat->GetBlockCount()];

      for (unsigned int ii = 0; ii < iGlobLatDat->GetBlockCount(); ++ii)
      {
        readBlock[ii] = false;
      }

      // Read a block in if it has fluid sites and is to live on the current processor. Also read
      // in any neighbours with fluid sites.
      for (unsigned int blockI = 0; blockI < iGlobLatDat->GetXBlockCount(); ++blockI)
      {
        for (unsigned int blockJ = 0; blockJ < iGlobLatDat->GetYBlockCount(); ++blockJ)
        {
          for (unsigned int blockK = 0; blockK < iGlobLatDat->GetZBlockCount(); ++blockK)
          {
            unsigned int lBlockId = iGlobLatDat->GetBlockIdFromBlockCoords(blockI, blockJ, blockK);

            if ( (bytesPerBlock[lBlockId] == 0) || (unitForEachBlock[lBlockId] != (int) localRank))
            {
              continue;
            }

            // Read in all neighbouring blocks with fluid sites.
            for (unsigned int neighI = util::NumericalFunctions::max<int>(0, blockI - 1); (neighI
                <= (blockI + 1)) && (neighI < iGlobLatDat->GetXBlockCount()); ++neighI)
            {
              for (unsigned int neighJ = util::NumericalFunctions::max<int>(0, blockJ - 1); (neighJ
                  <= (blockJ + 1)) && (neighJ < iGlobLatDat->GetYBlockCount()); ++neighJ)
              {
                for (unsigned int neighK = util::NumericalFunctions::max<int>(0, blockK - 1); (neighK
                    <= (blockK + 1)) && (neighK < iGlobLatDat->GetZBlockCount()); ++neighK)
                {
                  unsigned int lNeighId = iGlobLatDat->GetBlockIdFromBlockCoords(neighI, neighJ,
                                                                                 neighK);

                  if (bytesPerBlock[lNeighId] > 0)
                  {
                    readBlock[lNeighId] = true;
                  }
                }
              }
            }
          }
        }
      }

      // Each site can have at most
      // * an unsigned int (config)
      // * 4 doubles
      // * 14 further doubles (in D3Q15)
      unsigned int maxBytesPerBlock = (iGlobLatDat->SitesPerBlockVolumeUnit) * (4 * 1 + 8 * (4
          + D3Q15::NUMVECTORS - 1));
      const unsigned int BlocksToReadInOneGo = 10;
      char* readBuffer = new char[maxBytesPerBlock * BlocksToReadInOneGo];

      // This makes sure we do the right number of read operations.
      for (unsigned int readNum = 0; readNum
          <= (iGlobLatDat->GetBlockCount() / BlocksToReadInOneGo); ++readNum)
      {
        const unsigned int upperLimitBlockNumber =
            util::NumericalFunctions::min<unsigned int>(iGlobLatDat->GetBlockCount(), (readNum + 1)
                * BlocksToReadInOneGo);
        const unsigned int lowerLimitBlockNumber = readNum * BlocksToReadInOneGo;

        unsigned int bytesToRead = 0;
        for (unsigned int ii = lowerLimitBlockNumber; ii < upperLimitBlockNumber; ++ii)
        {
          bytesToRead += bytesPerBlock[ii];
        }

        MPI_Status lStatus;

        MPI_File_read_all(iFile, readBuffer, (int) bytesToRead, MPI_BYTE, &lStatus);

        io::XdrMemReader lReader(readBuffer, bytesToRead);

        for (lb::BlockCounter lBlock(iGlobLatDat, lowerLimitBlockNumber); lBlock
            < upperLimitBlockNumber; lBlock++)
        {
          if (readBlock[lBlock] && bytesPerBlock[lBlock] > 0)
          {
            iGlobLatDat->Blocks[lBlock].site_data
                = new unsigned int[iGlobLatDat->SitesPerBlockVolumeUnit];
            iGlobLatDat->Blocks[lBlock].ProcessorRankForEachBlockSite
                = new int[iGlobLatDat->SitesPerBlockVolumeUnit];

            int m = -1;

            for (unsigned int ii = 0; ii < iGlobLatDat->GetBlockSize(); ii++)
            {
              unsigned int site_i = lBlock.GetICoord(ii);

              for (unsigned int jj = 0; jj < iGlobLatDat->GetBlockSize(); jj++)
              {
                unsigned int site_j = lBlock.GetJCoord(jj);

                for (unsigned int kk = 0; kk < iGlobLatDat->GetBlockSize(); kk++)
                {
                  unsigned int site_k = lBlock.GetKCoord(kk);

                  ++m;

                  unsigned int *site_type = &iGlobLatDat->Blocks[lBlock].site_data[m];
                  if (!lReader.readUnsignedInt(*site_type))
                  {
                    std::cout << "Error reading site type\n";
                  }

                  if ( (*site_type & SITE_TYPE_MASK) == hemelb::lb::SOLID_TYPE)
                  {
                    iGlobLatDat->Blocks[lBlock].ProcessorRankForEachBlockSite[m] = 1 << 30;
                    continue;
                  }
                  iGlobLatDat->Blocks[lBlock].ProcessorRankForEachBlockSite[m] = -1;

                  if (iGlobLatDat->GetCollisionType(*site_type) != FLUID)
                  {
                    // Neither solid nor simple fluid
                    if (iGlobLatDat->Blocks[lBlock].wall_data == NULL)
                    {
                      iGlobLatDat->Blocks[lBlock].wall_data
                          = new hemelb::lb::WallData[iGlobLatDat->SitesPerBlockVolumeUnit];
                    }

                    if (iGlobLatDat->GetCollisionType(*site_type) & INLET
                        || iGlobLatDat->GetCollisionType(*site_type) & OUTLET)
                    {
                      double temp;
                      // INLET or OUTLET or both.
                      // These values are the boundary normal and the boundary distance.
                      for (int l = 0; l < 3; l++)
                      {
                        if (!lReader.readDouble(temp))
                        {
                          std::cout << "Error reading boundary normals\n";
                        }
                      }

                      if (!lReader.readDouble(temp))
                      {
                        std::cout << "Error reading boundary distances\n";
                      }
                    }

                    if (iGlobLatDat->GetCollisionType(*site_type) & EDGE)
                    {
                      // EDGE bit set
                      for (int l = 0; l < 3; l++)
                      {
                        if (!lReader.readDouble(
                                                iGlobLatDat->Blocks[lBlock].wall_data[m].wall_nor[l]))
                        {
                          std::cout << "Error reading edge normal\n";
                        }
                      }

                      double temp;
                      if (!lReader.readDouble(temp))
                      {
                        std::cout << "Error reading edge distance\n";
                      }
                    }

                    for (unsigned int l = 0; l < (D3Q15::NUMVECTORS - 1); l++)
                    {
                      if (!lReader.readDouble(iGlobLatDat->Blocks[lBlock].wall_data[m].cut_dist[l]))
                      {
                        std::cout << "Error reading cut distances\n";
                      }
                    }
                  }
                } // kk
              } // jj
            } // ii
          }
          // Not being read.
          else
          {
            if (bytesPerBlock[lBlock] > 0)
            {
              if (!lReader.SetPosition(lReader.GetPosition() + bytesPerBlock[lBlock]))
              {
                std::cout << "Error setting stream position\n";
              }
            }
          }
        }
      }

      delete[] readBuffer;
      delete[] readBlock;
    }

    void TopologyReader::DivideBlocks(unsigned int currentUnit,
                                      unsigned int blocksPerUnit,
                                      unsigned int unassignedBlocks,
                                      unsigned int totalBlockCount,
                                      unsigned int unitCount,
                                      unsigned int* blocksOnEachUnit,
                                      int* unitForEachBlock,
                                      const unsigned int* fluidSitesPerBlock,
                                      const lb::GlobalLatticeData* iGlobLatDat)
    {
      bool *blockAssigned = new bool[totalBlockCount];

      for (unsigned int ii = 0; ii < totalBlockCount; ++ii)
      {
        blockAssigned[ii] = false;
      }

      std::vector<BlockLocation> *lBlockLocationA = new std::vector<BlockLocation>;
      std::vector<BlockLocation> *lBlockLocationB = new std::vector<BlockLocation>;

      int lBlockNumber = -1;

      // Domain Decomposition.  Pick a site. Set it to the rank we are
      // looking at. Find its neighbours and put those on the same
      // rank, then find the next-nearest neighbours, etc. until we
      // have a completely joined region, or there are enough fluid
      // sites on the rank.  In the former case, start again at
      // another site. In the latter case, move on to the next rank.
      // Do this until all sites are assigned to a rank. There is a
      // high chance of of all sites on a rank being joined.

      unsigned int lBlocksOnCurrentProc = 0;

      // Iterate over all blocks.
      for (unsigned int lBlockCoordI = 0; lBlockCoordI < iGlobLatDat->GetXBlockCount(); lBlockCoordI++)
      {
        for (unsigned int lBlockCoordJ = 0; lBlockCoordJ < iGlobLatDat->GetYBlockCount(); lBlockCoordJ++)
        {
          for (unsigned int lBlockCoordK = 0; lBlockCoordK < iGlobLatDat->GetZBlockCount(); lBlockCoordK++)
          {
            // Block number is the number of the block we're currently on.
            lBlockNumber++;

            // If the array of proc rank for each site is NULL, we're on an all-solid block.
            // Alternatively, if this block has already been assigned, move on.
            if (fluidSitesPerBlock[lBlockNumber] == 0)
            {
              unitForEachBlock[lBlockNumber] = -1;
              continue;
            }
            else if (blockAssigned[lBlockNumber])
            {
              continue;
            }

            // Assign this block to the current unit.
            blockAssigned[lBlockNumber] = true;
            unitForEachBlock[lBlockNumber] = currentUnit;

            ++lBlocksOnCurrentProc;

            // Record the location of this initial site.
            lBlockLocationA->clear();
            BlockLocation lNew;
            lNew.i = lBlockCoordI;
            lNew.j = lBlockCoordJ;
            lNew.k = lBlockCoordK;
            lBlockLocationA->push_back(lNew);

            // The subdomain can grow.
            bool lIsRegionGrowing = true;

            // While the region can grow (i.e. it is not bounded by solids or visited
            // sites), and we need more sites on this particular rank.
            while (lBlocksOnCurrentProc < blocksPerUnit && lIsRegionGrowing)
            {
              lBlockLocationB->clear();

              // Sites added to the edge of the mClusters during the iteration.
              lIsRegionGrowing = false;

              // For sites on the edge of the domain (sites_a), deal with the neighbours.
              for (unsigned int index_a = 0; index_a < lBlockLocationA->size()
                  && lBlocksOnCurrentProc < blocksPerUnit; index_a++)
              {
                lNew = lBlockLocationA->operator [](index_a);

                for (unsigned int l = 1; l < D3Q15::NUMVECTORS && lBlocksOnCurrentProc
                    < blocksPerUnit; l++)
                {
                  // Record neighbour location.
                  int neigh_i = lNew.i + D3Q15::CX[l];
                  int neigh_j = lNew.j + D3Q15::CY[l];
                  int neigh_k = lNew.k + D3Q15::CZ[l];

                  // Move on if neighbour is outside the bounding box.
                  if (neigh_i == -1 || neigh_i == (int) iGlobLatDat->GetXBlockCount())
                    continue;
                  if (neigh_j == -1 || neigh_j == (int) iGlobLatDat->GetYBlockCount())
                    continue;
                  if (neigh_k == -1 || neigh_k == (int) iGlobLatDat->GetZBlockCount())
                    continue;

                  // Move on if the neighbour is in a block of solids (in which case
                  // the pointer to ProcessorRankForEachBlockSite is NULL) or it is solid or has already
                  // been assigned to a rank (in which case ProcessorRankForEachBlockSite != -1).  ProcessorRankForEachBlockSite
                  // was initialized in lbmReadConfig in io.cc.

                  unsigned int neighBlockId = iGlobLatDat->GetBlockIdFromBlockCoords(neigh_i,
                                                                                     neigh_j,
                                                                                     neigh_k);

                  // Don't use this block if it has no fluid sites, or if it has already been assigned to a processor.
                  if (fluidSitesPerBlock[neighBlockId] == 0 || blockAssigned[neighBlockId])
                  {
                    continue;
                  }

                  // Set the rank for a neighbour and update the fluid site counters.
                  blockAssigned[neighBlockId] = true;
                  unitForEachBlock[neighBlockId] = currentUnit;
                  lBlocksOnCurrentProc++;

                  // Neighbour was found, so the region can grow.
                  lIsRegionGrowing = true;

                  // Record the location of the neighbour.
                  BlockLocation lNewB;
                  lNewB.i = neigh_i;
                  lNewB.j = neigh_j;
                  lNewB.k = neigh_k;
                  lBlockLocationB->push_back(lNewB);
                }
              }
              // When the new layer of edge sites has been found, swap the buffers for
              // the current and new layers of edge sites.
              std::vector<BlockLocation> *tempP = lBlockLocationA;
              lBlockLocationA = lBlockLocationB;
              lBlockLocationB = tempP;
            }

            // If we have enough sites, we have finished.
            if (lBlocksOnCurrentProc >= blocksPerUnit)
            {
              blocksOnEachUnit[currentUnit] = lBlocksOnCurrentProc;

              ++currentUnit;

              unassignedBlocks -= lBlocksOnCurrentProc;
              blocksPerUnit = (int) ceil((double) unassignedBlocks / (double) (unitCount
                  - currentUnit));

              lBlocksOnCurrentProc = 0;
            }
            // If not, we have to start growing a different region for the same rank:
            // region expansions could get trapped.

          } // Block co-ord k
        } // Block co-ord j
      } // Block co-ord i

      delete lBlockLocationA;
      delete lBlockLocationB;

      delete[] blockAssigned;
    }

    void TopologyReader::OptimiseDomainDecomposition(const unsigned int* sitesPerBlock,
                                                     const int* procForEachBlock,
                                                     const unsigned int* blockCountForEachProc,
                                                     const topology::NetworkTopology* iNetTop,
                                                     SimConfig* bSimConfig,
                                                     lb::LbmParameters* bLbmParams,
                                                     lb::GlobalLatticeData* bGlobLatDat)
    {
      /*
       *  Get an array of the site count on each processor.
       */
      int* sitesPerProc = new int[iNetTop->GetProcessorCount()];
      for (unsigned int ii = 0; ii < iNetTop->GetProcessorCount(); ++ii)
      {
        sitesPerProc[ii] = 0;
      }

      for (unsigned int ii = 0; ii < bGlobLatDat->GetBlockCount(); ++ii)
      {
        sitesPerProc[procForEachBlock[ii]] += sitesPerBlock[ii];
      }

      /*
       *  Create vertex distribution array.
       */
      int* vertexDistribution = new int[iNetTop->GetProcessorCount() + 1];
      vertexDistribution[0] = 0;

      for (unsigned int ii = 0; ii < iNetTop->GetProcessorCount(); ++ii)
      {
        vertexDistribution[ii + 1] = vertexDistribution[ii] + sitesPerProc[ii];
      }

      /*
       *  Create the adjacency count and list for all the local vertices.
       */
      // First, it'll be useful to create an array of the first site index on each block.

      // This needs to have all the blocks with ascending id based on which processor they're on.
      int* firstSiteOnProc = new int[iNetTop->GetProcessorCount()];
      for (unsigned int ii = 0; ii < iNetTop->GetProcessorCount(); ++ii)
      {
        firstSiteOnProc[ii] = vertexDistribution[ii];
      }

      int* firstSiteIndexPerBlock = new int[bGlobLatDat->GetBlockCount()];
      for (unsigned int ii = 0; ii < bGlobLatDat->GetBlockCount(); ++ii)
      {
        int proc = procForEachBlock[ii];
        if (proc < 0)
        {
          firstSiteIndexPerBlock[ii] = -1;
        }
        else
        {
          firstSiteIndexPerBlock[ii] = firstSiteOnProc[proc];
          firstSiteOnProc[proc] += sitesPerBlock[ii];
        }
      }

      int localVertexCount = vertexDistribution[iNetTop->GetLocalRank() + 1]
          - vertexDistribution[iNetTop->GetLocalRank()];

      int* adjacenciesPerVertex = new int[localVertexCount + 1];
      adjacenciesPerVertex[0] = 0;

      std::vector<int> lAdjacencies;

      unsigned int lFluidVertex = 0;
      int n = -1;
      for (unsigned int i = 0; i < bGlobLatDat->GetXSiteCount(); i += bGlobLatDat->GetBlockSize())
      {
        for (unsigned int j = 0; j < bGlobLatDat->GetYSiteCount(); j += bGlobLatDat->GetBlockSize())
        {
          for (unsigned int k = 0; k < bGlobLatDat->GetZSiteCount(); k
              += bGlobLatDat->GetBlockSize())
          {
            ++n;

            if (procForEachBlock[n] != (int) iNetTop->GetLocalRank())
            {
              continue;
            }

            hemelb::lb::BlockData *map_block_p = &bGlobLatDat->Blocks[n];

            if (map_block_p->site_data == NULL)
            {
              continue;
            }

            int m = -1;

            // Iterate over sites within the block.
            for (unsigned int site_i = i; site_i < i + bGlobLatDat->GetBlockSize(); site_i++)
            {
              for (unsigned int site_j = j; site_j < j + bGlobLatDat->GetBlockSize(); site_j++)
              {
                for (unsigned int site_k = k; site_k < k + bGlobLatDat->GetBlockSize(); site_k++)
                {
                  ++m;

                  // Get site data, which is the number of the fluid site on this proc..
                  unsigned int site_map = map_block_p->site_data[m];

                  // Continue if it's a solid
                  if (map_block_p->ProcessorRankForEachBlockSite[m] == 1U << 30)
                  {
                    continue;
                  }

                  for (unsigned int l = 1; l < D3Q15::NUMVECTORS; l++)
                  {
                    // Work out positions of neighbours.
                    int neigh_i = site_i + D3Q15::CX[l];
                    int neigh_j = site_j + D3Q15::CY[l];
                    int neigh_k = site_k + D3Q15::CZ[l];

                    // Get the id of the processor which the neighbouring site lies on.
                    int *proc_id_p = bGlobLatDat->GetProcIdFromGlobalCoords(neigh_i, neigh_j,
                                                                            neigh_k);

                    if (proc_id_p == NULL || *proc_id_p == BIG_NUMBER2)
                    {
                      continue;
                    }

                    // We now do some faffery to find out the global fluid site id of this point
                    unsigned int neighBlockI = neigh_i >> bGlobLatDat->Log2BlockSize;
                    unsigned int neighBlockJ = neigh_j >> bGlobLatDat->Log2BlockSize;
                    unsigned int neighBlockK = neigh_k >> bGlobLatDat->Log2BlockSize;

                    unsigned int neighLocalSiteI = neigh_i - (neighBlockI
                        << bGlobLatDat->Log2BlockSize);
                    unsigned int neighLocalSiteJ = neigh_j - (neighBlockJ
                        << bGlobLatDat->Log2BlockSize);
                    unsigned int neighLocalSiteK = neigh_k - (neighBlockK
                        << bGlobLatDat->Log2BlockSize);

                    unsigned int neighBlockId = bGlobLatDat->GetBlockIdFromBlockCoords(neighBlockI,
                                                                                       neighBlockJ,
                                                                                       neighBlockK);

                    unsigned int neighGlobalSiteId = firstSiteIndexPerBlock[neighBlockId];

                    unsigned int localSiteId = ( ( (neighLocalSiteI << bGlobLatDat->Log2BlockSize)
                        + neighLocalSiteJ) << bGlobLatDat->Log2BlockSize) + neighLocalSiteK;

                    for (unsigned int neighSite = 0; neighSite
                        < bGlobLatDat->SitesPerBlockVolumeUnit; ++neighSite)
                    {
                      if (neighSite == localSiteId)
                      {
                        break;
                      }
                      else if (bGlobLatDat->Blocks[neighBlockId].ProcessorRankForEachBlockSite[neighSite]
                          != 1U << 30)
                      {
                        ++neighGlobalSiteId;
                      }
                    }

                    lAdjacencies.push_back(neighGlobalSiteId);
                  }

                  adjacenciesPerVertex[lFluidVertex + 1] = lAdjacencies.size();
                  ++lFluidVertex;
                }
              }
            }
          }
        }
      }

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
      /*
       unsigned long * lCumulativeSitesPerProc = new unsigned long[mTopologySize + 1];
       lCumulativeSitesPerProc[0] = 0;
       for (int ii = 0; ii < mTopologySize; ii++)
       {
       lCumulativeSitesPerProc[ii + 1] = lCumulativeSitesPerProc[ii] + iNumberSitesPerProc[ii];
       }

       //TODO

       int* siteCountPerProc = iNetTop->FluidSitesOnEachProcessor;
       int* cumulativeAdjsThisNode = new int[1
       + iNetTop->FluidSitesOnEachProcessor[iNetTop->GetLocalRank()]];

       cumulativeAdjsThisNode[0] = 0;
       std::vector<int> adjacencies;

       for (int ii = 0; ii < iNetTop->FluidSitesOnEachProcessor[iNetTop->GetLocalRank()]; ++ii)
       {

       }
       */
      // Required: zero-indexed site id for first site on each block.

      int weightFlag = 0;
      int numberingFlag = 0;
      int edgesCut = 0;
      int* partitionVector = new int[localVertexCount];
      int options[4] = { 0, 0, 0, 0 };

      MPI_Comm lComms = MPI_COMM_WORLD;

      unsigned int edgesCutBefore = 0;
      int myLowest = vertexDistribution[iNetTop->GetLocalRank()];
      int myHighest = vertexDistribution[iNetTop->GetLocalRank() + 1] - 1;

      for (unsigned int ii = 0; ii < lAdjacencies.size(); ++ii)
      {
        if (lAdjacencies[ii] < myLowest || lAdjacencies[ii] > myHighest)
        {
          ++edgesCutBefore;
        }
      }

      int desiredPartitionSize = iNetTop->GetProcessorCount() - 1;

      //      ParMETIS_RefineKway(vertexDistribution, adjacenciesPerVertex, &lAdjacencies[0], NULL, NULL,
      //                          &weightFlag, &numberingFlag, options, &edgesCut, partitionVector, &lComms);

      ParMETIS_PartKway(vertexDistribution, adjacenciesPerVertex, &lAdjacencies[0], NULL, NULL,
                        &weightFlag, &numberingFlag, &desiredPartitionSize, options, &edgesCut,
                        partitionVector, &lComms);

      // Because of the way Parmetis works, it will have assigned everything for the last processor to proc 0.
      for (int ii = 0; ii <= (myHighest - myLowest); ++ii)
      {
        ++ (partitionVector[ii]);
      }

      // Right. Let's count how many sites we're going to have to move.
      int moves = 0;
      std::vector<int> moveData;

      for (int ii = 0; ii <= (myHighest - myLowest); ++ii)
      {
        if (partitionVector[ii] != (int) iNetTop->GetLocalRank())
        {
          ++moves;
          moveData.push_back(myLowest + ii);
          moveData.push_back(partitionVector[ii]);
        }
      }

      // Spread this data around, so all processes now how many moves each process is doing.
      int* allMoves = new int[iNetTop->GetProcessorCount()];

      MPI_Allgather(&moves, 1, MPI_INT, allMoves, 1, MPI_INT, MPI_COMM_WORLD);

      // Count the total moves.
      int totalMoves = 0;

      for (unsigned int ii = 0; ii < iNetTop->GetProcessorCount(); ++ii)
      {
        totalMoves += allMoves[ii];
      }

      // Now share all the lists of moves.
      MPI_Datatype lMoveType;
      MPI_Type_contiguous(2, MPI_INT, &lMoveType);
      MPI_Type_commit(&lMoveType);

      int* movesList = new int[2 * totalMoves];
      int* offsets = new int[iNetTop->GetProcessorCount()];

      offsets[0] = 0;

      for (unsigned int ii = 1; ii < iNetTop->GetProcessorCount(); ++ii)
      {
        offsets[ii] = offsets[ii - 1] + allMoves[ii - 1];
      }

      MPI_Allgatherv(&moveData[0], moves, lMoveType, movesList, allMoves, offsets, lMoveType,
                     MPI_COMM_WORLD);

      int* newProcForEachBlock = new int[bGlobLatDat->GetBlockCount()];

      for (unsigned int lBlockNumber = 0; lBlockNumber < bGlobLatDat->GetBlockCount(); ++lBlockNumber)
      {
        newProcForEachBlock[lBlockNumber] = procForEachBlock[lBlockNumber];
      }

      // Hackily, set the proc for each block to be this proc whenever we have a proc on that block under the new scheme.
      unsigned int moveIndex = 0;

      for (unsigned int lFromProc = 0; lFromProc < iNetTop->GetProcessorCount(); ++lFromProc)
      {
        for (int lMoveNumber = 0; lMoveNumber < allMoves[lFromProc]; ++lMoveNumber)
        {
          int fluidSite = movesList[2 * moveIndex];
          int toProc = movesList[2 * moveIndex + 1];
          ++moveIndex;

          if (toProc == (int) iNetTop->GetLocalRank())
          {
            unsigned int fluidSiteBlock = 0;

            while ( (procForEachBlock[fluidSiteBlock] < 0)
                || (firstSiteIndexPerBlock[fluidSiteBlock] > fluidSite)
                || ( (firstSiteIndexPerBlock[fluidSiteBlock] + (int) sitesPerBlock[fluidSiteBlock])
                    <= fluidSite))
            {
              fluidSiteBlock++;
            }

            newProcForEachBlock[fluidSiteBlock] = (int) iNetTop->GetLocalRank();
          }
        }
      }

      // TODO USE A COMMUNICATOR JUST FOR THIS CLASS.
      MPI_Type_free(&lMoveType);

      MPI_File lFile;
      MPI_Info lFileInfo;
      int lError;

      // Open the file using the MPI parallel I/O interface at the path
      // given, in read-only mode.
      MPI_Info_create(&lFileInfo);

      // Create hints about how we'll read the file. See Chapter 13, page 400 of the MPI 2.2 spec.
      //
      MPI_Info_set(lFileInfo, "access_style", "read_once,sequential");
      MPI_Info_set(lFileInfo, "collective_buffering", "true");

      lError = MPI_File_open(MPI_COMM_WORLD, &bSimConfig->DataFilePath[0], MPI_MODE_RDONLY,
                             lFileInfo, &lFile);

      if (lError != 0)
      {
        fprintf(stderr, "Unable to open file %s [rank %i], exiting\n",
                bSimConfig->DataFilePath.c_str(), mGlobalRank);
        fflush(0x0);
        exit(0x0);
      }

      ReadPreamble(lFile, bLbmParams, bGlobLatDat);

      unsigned int* bytesPerBlock = new unsigned int[bGlobLatDat->GetBlockCount()];
      unsigned int* newSitePerBlock = new unsigned int[bGlobLatDat->GetBlockCount()];

      ReadHeader(lFile, bGlobLatDat->GetBlockCount(), newSitePerBlock, bytesPerBlock);

      ReadInLocalBlocks(lFile, bytesPerBlock, newProcForEachBlock, iNetTop->GetLocalRank(),
                        bGlobLatDat);

      // Set the proc rank to what it originally was.
      for (unsigned int lBlock = 0; lBlock < bGlobLatDat->GetBlockCount(); ++lBlock)
      {
        int lOriginalProc = procForEachBlock[lBlock];

        if (bGlobLatDat->Blocks[lBlock].ProcessorRankForEachBlockSite != NULL)
        {
          for (unsigned int lSiteIndex = 0; lSiteIndex < bGlobLatDat->SitesPerBlockVolumeUnit; ++lSiteIndex)
          {
            if (bGlobLatDat->Blocks[lBlock].ProcessorRankForEachBlockSite[lSiteIndex] != (1U << 30))
            {
              // This should be what it originally was...
              bGlobLatDat->Blocks[lBlock].ProcessorRankForEachBlockSite[lSiteIndex] = lOriginalProc;
            }
          }
        }
      }

      // Then go through and implement the ParMetis partition.
      moveIndex = 0;

      for (unsigned int lFromProc = 0; lFromProc < iNetTop->GetProcessorCount(); ++lFromProc)
      {
        for (int lMoveNumber = 0; lMoveNumber < allMoves[lFromProc]; ++lMoveNumber)
        {
          int fluidSite = movesList[2 * moveIndex];
          int toProc = movesList[2 * moveIndex + 1];
          ++moveIndex;

          unsigned int fluidSiteBlock = 0;

          while ( (procForEachBlock[fluidSiteBlock] < 0) || (firstSiteIndexPerBlock[fluidSiteBlock]
              > fluidSite) || ( (firstSiteIndexPerBlock[fluidSiteBlock]
              + (int) sitesPerBlock[fluidSiteBlock]) <= fluidSite))
          {
            fluidSiteBlock++;
          }

          if (bGlobLatDat->Blocks[fluidSiteBlock].ProcessorRankForEachBlockSite != NULL)
          {
            int fluidSitesToPass = fluidSite - firstSiteIndexPerBlock[fluidSiteBlock];

            unsigned int lSiteIndex = 0;

            while (true)
            {
              if (bGlobLatDat->Blocks[fluidSiteBlock].ProcessorRankForEachBlockSite[lSiteIndex]
                  != (1U << 30))
              {
                fluidSitesToPass--;
              }
              if (fluidSitesToPass < 0)
              {
                break;
              }
              lSiteIndex++;
            }

            bGlobLatDat->Blocks[fluidSiteBlock].ProcessorRankForEachBlockSite[lSiteIndex] = toProc;
          }
        }
      }

      MPI_File_close(&lFile);
      MPI_Info_free(&lFileInfo);

      delete[] bytesPerBlock;
      delete[] newSitePerBlock;
      delete[] newProcForEachBlock;
    }

  }
}
