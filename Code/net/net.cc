/*! \file net.cc
 \brief In this file the functions useful to discover the topology used and
 to create and delete the domain decomposition and the various
 buffers are defined.
 */

#include <cstdlib>
#include <cmath>
#include <cstdio>

#include "net/net.h"
#include "util/utilityFunctions.h"

namespace hemelb
{
  namespace net
  {

    void Net::EnsureEnoughRequests(unsigned int count)
    {
      if (mRequests.size() < count)
      {
        int deficit = count - mRequests.size();
        for (int ii = 0; ii < deficit; ii++)
        {
          mRequests.push_back(MPI_Request());
          mStatuses.push_back(MPI_Status());
        }
      }
    }

    /*!
     This is called from the main function.  First function to deal with processors.
     The domain partitioning technique and the management of the
     buffers useful for the inter-processor communications are
     implemented in this function.  The domain decomposition is based
     on a graph growing partitioning technique.
     */
    int* Net::Initialise(hemelb::lb::GlobalLatticeData &iGlobLatDat,
                         hemelb::lb::LocalLatticeData* &bLocalLatDat)
    {
      // Create a map between the two-level data representation and the 1D
      // compact one is created here.

      // This rank's site data.
      unsigned int
          *lThisRankSiteData =
              new unsigned int[mNetworkTopology->FluidSitesOnEachProcessor[mNetworkTopology->GetLocalRank()]];

      GetThisRankSiteData(iGlobLatDat, lThisRankSiteData);

      // The numbers of inter- and intra-machine neighbouring processors,
      // interface-dependent and independent fluid sites and shared
      // distribution functions of the reference processor are calculated
      // here.  neigh_proc is a static array that is declared in config.h.

      CountCollisionTypes(bLocalLatDat, iGlobLatDat, lThisRankSiteData);

      // the precise interface-dependent data (interface-dependent fluid
      // site locations and identifiers of the distribution functions
      // streamed between different partitions) are collected and the
      // buffers needed for the communications are set from here

      short int *f_data = new short int[4 * mNetworkTopology->TotalSharedFs];

      // Allocate the index in which to put the distribution functions received from the other
      // process.
      int* f_recv_iv = new int[mNetworkTopology->TotalSharedFs];

      short int ** lSharedFLocationForEachProc =
          new short int*[mNetworkTopology->NeighbouringProcs.size()];

      // Reset to zero again.
      mNetworkTopology->TotalSharedFs = 0;

      // Set the remaining neighbouring processor data.
      for (unsigned int n = 0; n < mNetworkTopology->NeighbouringProcs.size(); n++)
      {
        // f_data compacted according to number of shared f_s on each process.
        // f_data will be set later.


        // Site co-ordinates for each of the shared distribution, and the number of
        // the corresponding direction vector. Array is 4 elements for each shared distribution.
        lSharedFLocationForEachProc[n] = &f_data[mNetworkTopology->TotalSharedFs << 2];

        // Pointing to a few things, but not setting any variables.
        // FirstSharedF points to start of shared_fs.
        mNetworkTopology->NeighbouringProcs[n]->FirstSharedF
            = bLocalLatDat->GetLocalFluidSiteCount() * D3Q15::NUMVECTORS + 1
                + mNetworkTopology->TotalSharedFs;

        mNetworkTopology->TotalSharedFs += mNetworkTopology->NeighbouringProcs[n]->SharedFCount;
      }

      mNetworkTopology->NeighbourIndexFromProcRank
          = new short int[mNetworkTopology->GetProcessorCount()];

      for (unsigned int m = 0; m < mNetworkTopology->GetProcessorCount(); m++)
      {
        mNetworkTopology->NeighbourIndexFromProcRank[m] = -1;
      }
      // Get neigh_proc_index from ProcessorRankForEachBlockSite.
      for (unsigned int m = 0; m < mNetworkTopology->NeighbouringProcs.size(); m++)
      {
        mNetworkTopology->NeighbourIndexFromProcRank[mNetworkTopology->NeighbouringProcs[m]->Rank]
            = m;
      }

      InitialiseNeighbourLookup(bLocalLatDat, lSharedFLocationForEachProc, lThisRankSiteData,
                                iGlobLatDat);

      delete[] lThisRankSiteData;

      InitialisePointToPointComms(lSharedFLocationForEachProc);

      int f_count = bLocalLatDat->GetLocalFluidSiteCount() * D3Q15::NUMVECTORS;

      int sharedSitesSeen = 0;

      for (unsigned int m = 0; m < mNetworkTopology->NeighbouringProcs.size(); m++)
      {
        hemelb::topology::NeighbouringProcessor *neigh_proc_p =
            mNetworkTopology->NeighbouringProcs[m];

        for (int n = 0; n < neigh_proc_p->SharedFCount; n++)
        {
          // Get coordinates and direction of the distribution function to be sent to another process.
          short int *f_data_p = &lSharedFLocationForEachProc[m][n * 4];
          short int i = f_data_p[0];
          short int j = f_data_p[1];
          short int k = f_data_p[2];
          short int l = f_data_p[3];

          // Correct so that each process has the correct coordinates.
          if (neigh_proc_p->Rank < mNetworkTopology->GetLocalRank())
          {
            i += D3Q15::CX[l];
            j += D3Q15::CY[l];
            k += D3Q15::CZ[l];
            l = D3Q15::INVERSEDIRECTIONS[l];
          }

          // Get the fluid site number of site that will send data to another process.
          unsigned int site_map = *iGlobLatDat.GetSiteData(i, j, k);

          // Set f_id to the element in the send buffer that we put the updated
          // distribution functions in.
          bLocalLatDat->SetNeighbourLocation(site_map, l, ++f_count);

          // Set the place where we put the received distribution functions, which is
          // f_new[number of fluid site that sends, inverse direction].
          f_recv_iv[sharedSitesSeen] = site_map * D3Q15::NUMVECTORS + D3Q15::INVERSEDIRECTIONS[l];
          ++sharedSitesSeen;
        }
      }
      // neigh_prc->f_data was only set as a pointer to f_data, not allocated.  In this line, we
      // are freeing both of those.
      delete[] f_data;

      // Delete the array in which we kept the shared f locations. Don't delete subarrays - these
      // are pointers to elsewhere.
      delete[] lSharedFLocationForEachProc;

      return f_recv_iv;
    }

    void Net::InitialisePointToPointComms(short int **& lSharedFLocationForEachProc)
    {
      EnsureEnoughRequests(mNetworkTopology->NeighbouringProcs.size());

      // point-to-point communications are performed to match data to be
      // sent to/receive from different partitions; in this way, the
      // communication of the locations of the interface-dependent fluid
      // sites and the identifiers of the distribution functions which
      // propagate to different partitions is avoided (only their values
      // will be communicated). It's here!
      // Allocate the request variable.
      for (unsigned int m = 0; m < mNetworkTopology->NeighbouringProcs.size(); m++)
      {
        hemelb::topology::NeighbouringProcessor *neigh_proc_p =
            mNetworkTopology->NeighbouringProcs[m];
        // One way send receive.  The lower numbered mNetworkTopology->ProcessorCount send and the higher numbered ones receive.
        // It seems that, for each pair of processors, the lower numbered one ends up with its own
        // edge sites and directions stored and the higher numbered one ends up with those on the
        // other processor.
        if (neigh_proc_p->Rank > mNetworkTopology->GetLocalRank())
        {
          MPI_Isend(&lSharedFLocationForEachProc[m][0], neigh_proc_p->SharedFCount * 4, MPI_SHORT,
                    neigh_proc_p->Rank, 10, MPI_COMM_WORLD, &mRequests[m]);
        }
        else
        {
          MPI_Irecv(&lSharedFLocationForEachProc[m][0], neigh_proc_p->SharedFCount * 4, MPI_SHORT,
                    neigh_proc_p->Rank, 10, MPI_COMM_WORLD, &mRequests[m]);
        }
      }

      MPI_Waitall(mNetworkTopology->NeighbouringProcs.size(), &mRequests[0], &mStatuses[0]);
    }

    void Net::GetThisRankSiteData(const hemelb::lb::GlobalLatticeData &iGlobLatDat,
                                  unsigned int* &bThisRankSiteData)
    {
      // Array of booleans to store whether any sites on a block are fluid
      // sites residing on this rank.
      bool *lBlockIsOnThisRank = new bool[iGlobLatDat.GetBlockCount()];
      // Initialise to false.
      for (unsigned int n = 0; n < iGlobLatDat.GetBlockCount(); n++)
      {
        lBlockIsOnThisRank[n] = false;
      }

      int lSiteIndexOnProc = 0;

      for (unsigned int lBlockNumber = 0; lBlockNumber < iGlobLatDat.GetBlockCount(); lBlockNumber++)
      {
        hemelb::lb::BlockData * lCurrentDataBlock = &iGlobLatDat.Blocks[lBlockNumber];

        // If we are in a block of solids, move to the next block.
        if (lCurrentDataBlock->site_data == NULL)
        {
          continue;
        }

        // If we have some fluid sites, point to mProcessorsForEachBlock and map_block.
        hemelb::lb::BlockData * proc_block_p = &iGlobLatDat.Blocks[lBlockNumber];

        // lCurrentDataBlock.site_data is set to the fluid site identifier on this rank or (1U << 31U) if a site is solid
        // or not on this rank.  site_data is indexed by fluid site identifier and set to the site_data.
        for (unsigned int lSiteIndexWithinBlock = 0; lSiteIndexWithinBlock
            < iGlobLatDat.SitesPerBlockVolumeUnit; lSiteIndexWithinBlock++)
        {
          if ((int) mNetworkTopology->GetLocalRank()
              == proc_block_p->ProcessorRankForEachBlockSite[lSiteIndexWithinBlock])
          {
            // If the current site is non-solid, copy the site data into the array for
            // this rank (in the whole-processor location), then set the site data
            // for this site within the current block to be the site index over the whole
            // processor.
            if ( (lCurrentDataBlock->site_data[lSiteIndexWithinBlock] & SITE_TYPE_MASK)
                != hemelb::lb::SOLID_TYPE)
            {
              bThisRankSiteData[lSiteIndexOnProc]
                  = lCurrentDataBlock->site_data[lSiteIndexWithinBlock];
              lCurrentDataBlock->site_data[lSiteIndexWithinBlock] = lSiteIndexOnProc;
              ++lSiteIndexOnProc;
            }
            else
            {
              // If this is a solid, set the site data on the current block to
              // some massive value.
              lCurrentDataBlock->site_data[lSiteIndexWithinBlock] = (1U << 31U);
            }
            // Set the array to notify that the current block has sites on this
            // rank.
            lBlockIsOnThisRank[lBlockNumber] = true;
          }
          // If this site is not on the current processor, set its whole processor
          // index within the per-block store to a nonsense value.
          else
          {
            lCurrentDataBlock->site_data[lSiteIndexWithinBlock] = (1U << 31U);
          }
        }
      }

      // If we are in a block of solids, we set map_block[n].site_data to NULL.
      for (unsigned int n = 0; n < iGlobLatDat.GetBlockCount(); n++)
      {
        if (lBlockIsOnThisRank[n])
        {
          continue;
        }

        if (iGlobLatDat.Blocks[n].site_data != NULL)
        {
          delete[] iGlobLatDat.Blocks[n].site_data;
          iGlobLatDat.Blocks[n].site_data = NULL;
        }

        if (iGlobLatDat.Blocks[n].wall_data != NULL)
        {
          delete[] iGlobLatDat.Blocks[n].wall_data;
          iGlobLatDat.Blocks[n].wall_data = NULL;
        }
      }
      delete[] lBlockIsOnThisRank;
    }

    void Net::CountCollisionTypes(hemelb::lb::LocalLatticeData * bLocalLatDat,
                                  const hemelb::lb::GlobalLatticeData & iGlobLatDat,
                                  const unsigned int * lThisRankSiteData)
    {
      // Initialise various things to 0.
      bLocalLatDat->my_inner_sites = 0;

      for (int m = 0; m < COLLISION_TYPES; m++)
      {
        bLocalLatDat->my_inter_collisions[m] = 0;
        bLocalLatDat->my_inner_collisions[m] = 0;
      }

      mNetworkTopology->TotalSharedFs = 0;

      int lSiteIndexOnProc = 0;

      int n = -1;

      // Iterate over all blocks in site units
      for (unsigned int i = 0; i < iGlobLatDat.GetXSiteCount(); i += iGlobLatDat.GetBlockSize())
      {
        for (unsigned int j = 0; j < iGlobLatDat.GetYSiteCount(); j += iGlobLatDat.GetBlockSize())
        {
          for (unsigned int k = 0; k < iGlobLatDat.GetZSiteCount(); k += iGlobLatDat.GetBlockSize())
          {
            hemelb::lb::BlockData * map_block_p = &iGlobLatDat.Blocks[++n];

            if (map_block_p->site_data == NULL)
            {
              continue;
            }

            int m = -1;

            // Iterate over all sites within the current block.
            for (unsigned int site_i = i; site_i < i + iGlobLatDat.GetBlockSize(); site_i++)
            {
              for (unsigned int site_j = j; site_j < j + iGlobLatDat.GetBlockSize(); site_j++)
              {
                for (unsigned int site_k = k; site_k < k + iGlobLatDat.GetBlockSize(); site_k++)
                {
                  m++;
                  // If the site is not on this processor, continue.
                  if ((int) mNetworkTopology->GetLocalRank()
                      != map_block_p->ProcessorRankForEachBlockSite[m])
                  {
                    continue;
                  }

                  bool lIsInnerSite = true;

                  // Iterate over all direction vectors.
                  for (unsigned int l = 1; l < D3Q15::NUMVECTORS; l++)
                  {
                    // Find the neighbour site co-ords in this direction.
                    int neigh_i = site_i + D3Q15::CX[l];
                    int neigh_j = site_j + D3Q15::CY[l];
                    int neigh_k = site_k + D3Q15::CZ[l];

                    // Find the processor Id for that neighbour.
                    int *proc_id_p = iGlobLatDat.GetProcIdFromGlobalCoords(neigh_i, neigh_j,
                                                                           neigh_k);

                    // Move on if the neighbour is in a block of solids (in which case
                    // the pointer to ProcessorRankForEachBlockSite is NULL) or it is solid (in which case ProcessorRankForEachBlockSite ==
                    // BIG_NUMBER2) or the neighbour is also on this rank.  ProcessorRankForEachBlockSite was initialized
                    // in lbmReadConfig in io.cc.
                    if (proc_id_p == NULL || (int) mNetworkTopology->GetLocalRank() == (*proc_id_p)
                        || *proc_id_p == (BIG_NUMBER2))
                    {
                      continue;
                    }

                    lIsInnerSite = false;

                    // The first time, net_neigh_procs = 0, so
                    // the loop is not executed.
                    bool flag = true;

                    // Iterate over neighbouring processors until we find the one with the
                    // neighbouring site on it.
                    int lNeighbouringProcs = mNetworkTopology->NeighbouringProcs.size();
                    for (int mm = 0; mm < lNeighbouringProcs && flag; mm++)
                    {
                      // Check whether the rank for a particular neighbour has already been
                      // used for this processor.  If it has, set flag to zero.
                      hemelb::topology::NeighbouringProcessor * neigh_proc_p =
                          mNetworkTopology->NeighbouringProcs[mm];

                      // If ProcessorRankForEachBlockSite is equal to a neigh_proc that has alredy been listed.
                      if (*proc_id_p == (int) neigh_proc_p->Rank)
                      {
                        flag = false;
                        ++neigh_proc_p->SharedFCount;
                        ++mNetworkTopology->TotalSharedFs;
                      }
                    }
                    // If flag is 1, we didn't find a neighbour-proc with the neighbour-site on it
                    // so we need a new neighbouring processor.
                    if (flag)
                    {
                      // Store rank of neighbour in >neigh_proc[neigh_procs]
                      hemelb::topology::NeighbouringProcessor * lNewNeighbour =
                          new hemelb::topology::NeighbouringProcessor();
                      lNewNeighbour->SharedFCount = 1;
                      lNewNeighbour->Rank = *proc_id_p;
                      mNetworkTopology->NeighbouringProcs.push_back(lNewNeighbour);
                      ++mNetworkTopology->TotalSharedFs;
                    }
                  }

                  // Set the collision type data. map_block site data is renumbered according to
                  // fluid site numbers within a particular collision type.

                  int l = -1;

                  switch (iGlobLatDat.GetCollisionType(lThisRankSiteData[lSiteIndexOnProc]))
                  {
                    case FLUID:
                      l = 0;
                      break;
                    case EDGE:
                      l = 1;
                      break;
                    case INLET:
                      l = 2;
                      break;
                    case OUTLET:
                      l = 3;
                      break;
                    case (INLET | EDGE):
                      l = 4;
                      break;
                    case (OUTLET | EDGE):
                      l = 5;
                      break;
                  }

                  ++lSiteIndexOnProc;

                  if (lIsInnerSite)
                  {
                    ++bLocalLatDat->my_inner_sites;

                    if (l == 0)
                    {
                      map_block_p->site_data[m] = bLocalLatDat->my_inner_collisions[l];
                    }
                    else
                    {
                      map_block_p->site_data[m] = 50000000 * (10 + (l - 1))
                          + bLocalLatDat->my_inner_collisions[l];
                    }
                    ++bLocalLatDat->my_inner_collisions[l];
                  }
                  else
                  {
                    if (l == 0)
                    {
                      map_block_p->site_data[m] = 1000000000 + bLocalLatDat->my_inter_collisions[l];
                    }
                    else
                    {
                      map_block_p->site_data[m] = 50000000 * (20 + l)
                          + bLocalLatDat->my_inter_collisions[l];
                    }
                    ++bLocalLatDat->my_inter_collisions[l];
                  }
                }
              }
            }
          }
        }
      }

      int collision_offset[2][COLLISION_TYPES];
      // Calculate the number of each type of collision.
      collision_offset[0][0] = 0;

      for (unsigned int l = 1; l < COLLISION_TYPES; l++)
      {
        collision_offset[0][l] = collision_offset[0][l - 1] + bLocalLatDat->my_inner_collisions[l
            - 1];
      }
      collision_offset[1][0] = bLocalLatDat->my_inner_sites;
      for (unsigned int l = 1; l < COLLISION_TYPES; l++)
      {
        collision_offset[1][l] = collision_offset[1][l - 1] + bLocalLatDat->my_inter_collisions[l
            - 1];
      }

      // Iterate over blocks
      for (unsigned int n = 0; n < iGlobLatDat.GetBlockCount(); n++)
      {
        hemelb::lb::BlockData *map_block_p = &iGlobLatDat.Blocks[n];

        // If we are in a block of solids, continue.
        if (map_block_p->site_data == NULL)
        {
          continue;
        }

        // Iterate over sites within the block.
        for (unsigned int m = 0; m < iGlobLatDat.SitesPerBlockVolumeUnit; m++)
        {
          unsigned int *site_data_p = &map_block_p->site_data[m];

          // If the site is solid, continue.
          if (*site_data_p & (1U << 31U))
          {
            continue;
          }

          // 0th collision type for inner sites, so don't do anything.
          if (*site_data_p < 500000000)
          {
            continue;
          }

          // Renumber the sites in map_block so that the numbers are compacted together.  We have
          // collision offset to tell us when one collision type ends and another starts.
          for (unsigned int l = 1; l < COLLISION_TYPES; l++)
          {
            if (*site_data_p >= 50000000 * (10 + (l - 1)) && *site_data_p < 50000000 * (10 + l))
            {
              *site_data_p += collision_offset[0][l] - 50000000 * (10 + (l - 1));
              break;
            }
          }
          for (unsigned int l = 0; l < COLLISION_TYPES; l++)
          {
            if (*site_data_p >= 50000000 * (20 + l) && *site_data_p < 50000000 * (20 + (l + 1)))
            {
              *site_data_p += collision_offset[1][l] - 50000000 * (20 + l);
              break;
            }
          }
        }
      }

      bLocalLatDat->SetSharedSiteCount(mNetworkTopology->TotalSharedFs);
    }

    void Net::InitialiseNeighbourLookup(hemelb::lb::LocalLatticeData* bLocalLatDat,
                                        short int ** bSharedFLocationForEachProc,
                                        const unsigned int * iSiteDataForThisRank,
                                        const hemelb::lb::GlobalLatticeData & iGlobLatDat)
    {
      int n = -1;
      int lSiteIndexOnProc = 0;
      int * lFluidSitesHandledForEachProc = new int[mNetworkTopology->GetProcessorCount()];

      for (unsigned int ii = 0; ii < mNetworkTopology->GetProcessorCount(); ii++)
      {
        lFluidSitesHandledForEachProc[ii] = 0;
      }

      // Iterate over blocks in global co-ords.
      for (unsigned int i = 0; i < iGlobLatDat.GetXSiteCount(); i += iGlobLatDat.GetBlockSize())
      {
        for (unsigned int j = 0; j < iGlobLatDat.GetYSiteCount(); j += iGlobLatDat.GetBlockSize())
        {
          for (unsigned int k = 0; k < iGlobLatDat.GetZSiteCount(); k += iGlobLatDat.GetBlockSize())
          {
            n++;
            hemelb::lb::BlockData *map_block_p = &iGlobLatDat.Blocks[n];

            if (map_block_p->site_data == NULL)
            {
              continue;
            }

            int m = -1;

            // Iterate over sites within the block.
            for (unsigned int site_i = i; site_i < i + iGlobLatDat.GetBlockSize(); site_i++)
            {
              for (unsigned int site_j = j; site_j < j + iGlobLatDat.GetBlockSize(); site_j++)
              {
                for (unsigned int site_k = k; site_k < k + iGlobLatDat.GetBlockSize(); site_k++)
                {
                  // If a site is not on this process, continue.
                  m++;

                  if ((int) mNetworkTopology->GetLocalRank()
                      != map_block_p->ProcessorRankForEachBlockSite[m])
                  {
                    continue;
                  }

                  // Get site data, which is the number of the fluid site on this proc..
                  unsigned int site_map = map_block_p->site_data[m];

                  // Set neighbour location for the distribution component at the centre of
                  // this site.
                  bLocalLatDat->SetNeighbourLocation(site_map, 0, site_map * D3Q15::NUMVECTORS + 0);

                  for (unsigned int l = 1; l < D3Q15::NUMVECTORS; l++)
                  {
                    // Work out positions of neighbours.
                    int neigh_i = site_i + D3Q15::CX[l];
                    int neigh_j = site_j + D3Q15::CY[l];
                    int neigh_k = site_k + D3Q15::CZ[l];

                    // Get the id of the processor which the neighbouring site lies on.
                    int *proc_id_p = iGlobLatDat.GetProcIdFromGlobalCoords(neigh_i, neigh_j,
                                                                           neigh_k);

                    if (proc_id_p == NULL || *proc_id_p == BIG_NUMBER2)
                    {
                      // initialize f_id to the rubbish site.
                      bLocalLatDat->SetNeighbourLocation(site_map, l,
                                                         bLocalLatDat->GetLocalFluidSiteCount()
                                                             * D3Q15::NUMVECTORS);
                      continue;
                    }
                    // If on the same proc, set f_id of the
                    // current site and direction to the
                    // site and direction that it sends to.
                    // If we check convergence, the data for
                    // each site is split into that for the
                    // current and previous cycles.
                    else if ((int) mNetworkTopology->GetLocalRank() == *proc_id_p)
                    {

                      // Pointer to the neighbour.
                      const unsigned int *site_data_p = iGlobLatDat.GetSiteData(neigh_i, neigh_j,
                                                                                neigh_k);

                      bLocalLatDat->SetNeighbourLocation(site_map, l, *site_data_p
                          * D3Q15::NUMVECTORS + l);

                      continue;
                    }
                    else
                    {
                      short int neigh_proc_index =
                          mNetworkTopology->NeighbourIndexFromProcRank[*proc_id_p];

                      // This stores some coordinates.  We
                      // still need to know the site number.
                      // neigh_proc[ n ].f_data is now
                      // set as well, since this points to
                      // f_data.  Every process has data for
                      // its neighbours which say which sites
                      // on this process are shared with the
                      // neighbour.
                      short int
                          *f_data_p =
                              &bSharedFLocationForEachProc[neigh_proc_index][lFluidSitesHandledForEachProc[neigh_proc_index]
                                  << 2];
                      f_data_p[0] = site_i;
                      f_data_p[1] = site_j;
                      f_data_p[2] = site_k;
                      f_data_p[3] = l;
                      ++lFluidSitesHandledForEachProc[neigh_proc_index];
                    }
                  }

                  // This is used in Calculate BC in IO.
                  bLocalLatDat->mSiteData[site_map] = iSiteDataForThisRank[lSiteIndexOnProc];

                  if (iGlobLatDat.GetCollisionType(bLocalLatDat->mSiteData[site_map]) & EDGE)
                  {
                    bLocalLatDat->SetWallNormal(site_map,
                                                iGlobLatDat.Blocks[n].wall_data[m].wall_nor);

                    bLocalLatDat->SetDistanceToWall(site_map,
                                                    iGlobLatDat.Blocks[n].wall_data[m].cut_dist);
                  }
                  else
                  {
                    double lBigDistance[3];
                    for (unsigned int ii = 0; ii < 3; ii++)
                      lBigDistance[ii] = BIG_NUMBER;
                    bLocalLatDat->SetWallNormal(site_map, lBigDistance);
                  }
                  ++lSiteIndexOnProc;
                }
              }
            }
          }
        }
      }

      delete[] lFluidSitesHandledForEachProc;
    }

    void Net::Receive()
    {
      // Make sure the MPI datatypes have been created.
      EnsurePreparedToSendReceive();

      int m = 0;

      for (std::map<int, ProcComms*>::iterator it = mReceiveProcessorComms.begin(); it
          != mReceiveProcessorComms.end(); ++it)
      {
        MPI_Irecv(it->second->PointerList.front(), 1, it->second->Type, it->first, 10,
                  MPI_COMM_WORLD, &mRequests[m]);
        ++m;
      }
    }

    void Net::Send()
    {
      // Make sure the datatypes have been created.
      EnsurePreparedToSendReceive();

      int m = 0;

      for (std::map<int, ProcComms*>::iterator it = mSendProcessorComms.begin(); it
          != mSendProcessorComms.end(); ++it)
      {
        MPI_Isend(it->second->PointerList.front(), 1, it->second->Type, it->first, 10,
                  MPI_COMM_WORLD, &mRequests[mReceiveProcessorComms.size() + m]);

        ++m;
      }
    }

    void Net::Wait(hemelb::lb::LocalLatticeData *bLocalLatDat)
    {
      MPI_Waitall(mSendProcessorComms.size() + mReceiveProcessorComms.size(), &mRequests[0],
                  &mStatuses[0]);

      sendReceivePrepped = false;

      for (std::map<int, ProcComms*>::iterator it = mReceiveProcessorComms.begin(); it
          != mReceiveProcessorComms.end(); it++)
      {
        MPI_Type_free(&it->second->Type);
      }
      mReceiveProcessorComms.clear();

      for (std::map<int, ProcComms*>::iterator it = mSendProcessorComms.begin(); it
          != mSendProcessorComms.end(); it++)
      {
        MPI_Type_free(&it->second->Type);
      }
      mSendProcessorComms.clear();
    }

    // Helper function to get the ProcessorCommunications object, and create it if it doesn't exist yet.
    Net::ProcComms* Net::GetProcComms(int iRank, bool iIsSend)
    {
      std::map<int, ProcComms*>* lMap = iIsSend
        ? &mSendProcessorComms
        : &mReceiveProcessorComms;

      std::map<int, ProcComms*>::iterator lValue = lMap->find(iRank);
      ProcComms *lComms;
      if (lValue == lMap->end())
      {
        lComms = new ProcComms();
        lMap->insert(std::pair<int, ProcComms*>(iRank, lComms));
      }
      else
      {
        lComms = (lValue ->second);
      }
      return lComms;
    }

    // Helper functions to add ints to the list.
    void Net::AddToList(int* iNew, int iLength, ProcComms *bMetaData)
    {
      bMetaData->PointerList.push_back(iNew);
      bMetaData->LengthList.push_back(iLength);
      bMetaData->TypeList.push_back(MPI_INT);
    }

    // Helper functions to add doubles to the list.
    void Net::AddToList(double* iNew, int iLength, ProcComms *bMetaData)
    {
      bMetaData->PointerList.push_back(iNew);
      bMetaData->LengthList.push_back(iLength);
      bMetaData->TypeList.push_back(MPI_DOUBLE);
    }

    // Helper functions to add floats to the list.
    void Net::AddToList(float* iNew, int iLength, ProcComms *bMetaData)
    {
      bMetaData->PointerList.push_back(iNew);
      bMetaData->LengthList.push_back(iLength);
      bMetaData->TypeList.push_back(MPI_FLOAT);
    }

    // Makes sure the MPI_Datatypes for sending and receiving have been created for every neighbour.
    void Net::EnsurePreparedToSendReceive()
    {
      if (sendReceivePrepped)
      {
        return;
      }

      for (std::map<int, ProcComms*>::iterator it = mSendProcessorComms.begin(); it
          != mSendProcessorComms.end(); it++)
      {
        ProcComms* lThisPC = (*it).second;

        CreateMPIType(lThisPC);
      }

      for (std::map<int, ProcComms*>::iterator it = mReceiveProcessorComms.begin(); it
          != mReceiveProcessorComms.end(); it++)
      {
        ProcComms* lThisPC = (*it).second;

        CreateMPIType(lThisPC);
      }

      EnsureEnoughRequests(mReceiveProcessorComms.size() + mSendProcessorComms.size());

      sendReceivePrepped = true;
    }

    // Helper function to create a MPI derived datatype given a list of pointers, types and lengths.
    void Net::CreateMPIType(ProcComms *iMetaData)
    {
      MPI_Aint* displacements = new MPI_Aint[iMetaData->PointerList.size()];
      int* lengths = new int[iMetaData->PointerList.size()];
      MPI_Datatype* types = new MPI_Datatype[iMetaData->PointerList.size()];

      int lLocation = 0;

      for (std::vector<void*>::const_iterator it = iMetaData->PointerList.begin(); it
          != iMetaData->PointerList.end(); it++)
      {
        MPI_Get_address(*it, &displacements[lLocation]);
        ++lLocation;
      }

      for (int ii = iMetaData->PointerList.size() - 1; ii >= 0; ii--)
      {
        displacements[ii] -= displacements[0];

        lengths[ii] = iMetaData->LengthList[ii];
        types[ii] = iMetaData->TypeList[ii];
      }

      // Create the type and commit it.
      MPI_Type_create_struct(iMetaData->PointerList.size(), lengths, displacements, types,
                             &iMetaData->Type);
      MPI_Type_commit(&iMetaData->Type);

      delete[] displacements;
      delete[] lengths;
      delete[] types;
    }

    Net::Net(hemelb::topology::NetworkTopology * iNetworkTopology)
    {
      mNetworkTopology = iNetworkTopology;
      sendReceivePrepped = false;
    }

    /*!
     Free the allocated data.
     */
    Net::~Net()
    {
      if (sendReceivePrepped)
      {
        for (std::map<int, ProcComms*>::iterator it = mSendProcessorComms.begin(); it
            != mSendProcessorComms.end(); it++)
        {
          MPI_Type_free(&it->second->Type);
          delete it->second;
        }

        for (std::map<int, ProcComms*>::iterator it = mReceiveProcessorComms.begin(); it
            != mReceiveProcessorComms.end(); it++)
        {
          MPI_Type_free(&it->second->Type);
          delete it->second;
        }

      }
    }

  }
}
