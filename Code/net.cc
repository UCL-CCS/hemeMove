/*! \file net.cc
 \brief In this file the functions useful to discover the topology used and
 to create and delete the domain decomposition and the various
 buffers are defined.
 */

#include "lb.h"
#include "net.h"
#include "utilityFunctions.h"
#include <cstdlib>
#include <cmath>
#include <cstdio>

double Net::GetCutDistance(int iSiteIndex, int iDirection) const
{
  return cut_distances[iSiteIndex * (D3Q15::NUMVECTORS - 1) + iDirection - 1];
}

bool Net::HasBoundary(int iSiteIndex, int iDirection) const
{
  unsigned int lBoundaryConfig = (net_site_data[iSiteIndex]
      & BOUNDARY_CONFIG_MASK) >> BOUNDARY_CONFIG_SHIFT;
  return (lBoundaryConfig & (1U << (iDirection - 1))) != 0;
}

int Net::GetBoundaryId(int iSiteIndex) const
{
  return (net_site_data[iSiteIndex] & BOUNDARY_ID_MASK) >> BOUNDARY_ID_SHIFT;
}

double *Net::GetNormalToWall(int iSiteIndex) const
{
  return &net_site_nor[iSiteIndex * 3];
}

/*!
 Low level function that finds the pointer to the rank on which a
 particular site resides.  ProcessorRankForEachBlockSite is the only member of mProcessorsForEachBlock
 (member of Net) for the site at global coordinate (site_i, site_j,
 site_k).  If the site is in an empty block, return NULL.
 */
int * Net::GetProcIdFromGlobalCoords(int iSiteI, int iSiteJ, int iSiteK)
{
  // If the given site location is outside the bounding box return a NULL
  // pointer.
  if (iSiteI < 0 || iSiteI >= sites_x || iSiteJ < 0 || iSiteJ >= sites_y
      || iSiteK < 0 || iSiteK >= sites_z)
  {
    return NULL;
  }

  int ii, jj, kk; // Coordinates of a site within the block
  ProcBlock *proc_block_p; // Pointer to the block


  // Block identifiers (i, j, k) of the site (site_i, site_j, site_k)
  int i = iSiteI >> shift;
  int j = iSiteJ >> shift;
  int k = iSiteK >> shift;

  proc_block_p = &mProcessorsForEachBlock[ (i * blocks_y + j) * blocks_z + k];

  if (proc_block_p->ProcessorRankForEachBlockSite == NULL) // If an empty (solid) block is addressed
    return NULL;
  else
  {
    // Find site coordinates within the block
    ii = iSiteI - (i << shift);
    jj = iSiteJ - (j << shift);
    kk = iSiteK - (k << shift);

    // Return pointer to ProcessorRankForEachBlockSite[site] (the only member of
    // mProcessorsForEachBlock)
    return &proc_block_p->ProcessorRankForEachBlockSite[ ( ( (ii << shift) + jj)
        << shift) + kk];
  }
}

/*!
 Low level function that finds a pointer to site_data (the only
 member of map_block (member of net)) for the site eith global
 coordinate (site_i, site_j, site_k).  If the site is in an empty
 block, return NULL.
 */
unsigned int *Net::netSiteMapPointer(int site_i, int site_j, int site_k)
{
  int i, j, k; // Coordinates of a block
  int ii, jj, kk; // Coordinates of a site within the block
  DataBlock *map_block_p; // Pointer to the block

  if (site_i < 0 || site_i >= sites_x || site_j < 0 || site_j >= sites_y
      || site_k < 0 || site_k >= sites_z) // If site is out of the bounding box.
    return NULL;

  // Block identifiers (i, j, k) of the site (site_i, site_j, site_k)
  i = site_i >> shift;
  j = site_j >> shift;
  k = site_k >> shift;

  map_block_p = &map_block[ (i * blocks_y + j) * blocks_z + k];

  if (map_block_p->site_data == NULL) // if an empty (solid) block is addressed
    return NULL;
  else
  {
    // Find site coordinates within the block
    ii = site_i - (i << shift);
    jj = site_j - (j << shift);
    kk = site_k - (k << shift);

    // Return pointer to site_data[site]
    return &map_block_p->site_data[ ( ( (ii << shift) + jj) << shift) + kk];
  }
}

bool Net::IsCurrentProcRank(int iRank)
{
  return (mRank == iRank);
}

bool Net::IsCurrentProcTheIOProc()
{
  return mRank == 0;
}

char *Net::GetCurrentProcIdentifier()
{
  // Note that this has enough space for an integer of 9 digits, so we're good up to a billion
  // processors.
  char* lRet = new char[15];
  sprintf(lRet, "rank %i", mRank);
  return lRet;
}

//#undef MPICHX_TOPOLOGY_DEPTHS
#ifdef MPICHX_TOPOLOGY_DEPTHS

/*!
 If one has more than one machine. The topology discovery mechanism is implemented in this function
 */
int Net::netFindTopology (int *depths)
{

  int *depth, **color;
  int machine_id, flag, is_found;
  int i, j, sum;

  *depths = 0;

  err = MPI_Attr_get (MPI_COMM_WORLD, MPICHX_TOPOLOGY_DEPTHS, &depth, &flag);

  if (err != MPI_SUCCESS || flag == 0)
  {
    return 0;
  }

  err = MPI_Attr_get (MPI_COMM_WORLD, MPICHX_TOPOLOGY_COLORS, &color, &flag);

  if (err != MPI_SUCCESS || flag == 0)
  {
    return 0;
  }

  net_machines = 0;

  machine_id = new int [mProcessorCount];
  procs_per_machine = new int [mProcessorCount];

  for (i = 0; i < mProcessorCount; i++)
  {
    procs_per_machine[ i ] = 0;
  }
  for (i = 0; i < mProcessorCount; i++)
  {
    if (depth[ i ] != 4) continue;

    *depths = max(*depths, depth[ i ]);

    for (j = 0, is_found = 0; j < net_machines && is_found == 0; j++)
    {
      if (color[ i ][ 3 ] == machine_id[ j ])
      {
        is_found = 1;
        ++procs_per_machine[ machine_id[j] ];
      }
    }
    if (is_found == 1) continue;

    machine_id[ net_machines ] = color[ i ][ 3 ];
    ++procs_per_machine[ net_machines ];
    ++net_machines;
  }
  net_machines = max(1, net_machines);

  if (net_machines == 1)
  {
    for (i = 0; i < mProcessorCount; i++)
    {
      machine_id[ i ] = 0;
    }
    procs_per_machine[ 0 ] = mProcessorCount;
  }
  else
  {
    for (i = 0; i < mProcessorCount; i++)
    {
      sum = 0;
      machine_id = 0;

      is_found = 0;

      while (!is_found)
      {
        if (sum + procs_per_machine[ machine_id ] > i)
        {
          is_found = 1;
          continue;
        }
        sum += procs_per_machine[ machine_id ];
        ++machine_id;
      }
      machine_id[ i ] = machine_id;
    }
  }
  return 1;
}

#else

/*!
 If one has more than one machine. The topology discovery mechanism is implemented in this function
 */
int Net::netFindTopology(int *depths)
{
  // the machine is assumed to be only one if this function is
  // used instead of the previous one

  *depths = 1;

  net_machines = 1;

  machine_id = new int[mProcessorCount];
  procs_per_machine = new int[net_machines];

  for (int i = 0; i < mProcessorCount; i++)
  {
    machine_id[i] = 0;
  }
  procs_per_machine[0] = mProcessorCount;

  return 1;
}
#endif

// Returns the type of collision/streaming update for the fluid site
// with data "site_data".
unsigned int Net::GetCollisionType(unsigned int site_data)
{
  unsigned int boundary_type;

  if (site_data == FLUID_TYPE)
  {
    return FLUID;
  }
  boundary_type = site_data & SITE_TYPE_MASK;

  if (boundary_type == FLUID_TYPE)
  {
    return EDGE;
  }
  if (! (site_data & PRESSURE_EDGE_MASK))
  {
    if (boundary_type == INLET_TYPE)
    {
      return INLET;
    }
    else
    {
      return OUTLET;
    }
  }
  else
  {
    if (boundary_type == INLET_TYPE)
    {
      return INLET | EDGE;
    }
    else
    {
      return OUTLET | EDGE;
    }
  }
}

void Net::Abort()
{
#ifndef NOMPI
  err = MPI_Abort(MPI_COMM_WORLD, 1);
#else
  exit(1);
#endif
}

void Net::ThisNeedsRenaming(int & proc_count,
                            int & fluid_sites_per_unit,
                            int & unvisited_fluid_sites,
                            int marker,
                            int unitLevel)
{

  int sites_buffer_size = 10000;
  SiteLocation *site_location_a = new SiteLocation[sites_buffer_size];
  SiteLocation *site_location_b = new SiteLocation[sites_buffer_size];

  int lBlockNumber = -1;

  // Domain Decomposition.  Pick a site. Set it to the rank we are
  // looking at. Find its neighbours and put those on the same
  // rank, then find the next-nearest neighbours, etc. until we
  // have a completely joined region, or there are enough fluid
  // sites on the rank.  In the former case, start again at
  // another site. In the latter case, move on to the next rank.
  // Do this until all sites are assigned to a rank. There is a
  // high chance of of all sites on a rank being joined.

  // Iterate over all blocks.
  for (int lBlockCoordI = 0; lBlockCoordI < blocks_x; lBlockCoordI++)
  {
    for (int lBlockCoordJ = 0; lBlockCoordJ < blocks_y; lBlockCoordJ++)
    {
      for (int lBlockCoordK = 0; lBlockCoordK < blocks_z; lBlockCoordK++)
      {
        lBlockNumber++;

        // Point to a block of ProcessorRankForEachBlockSite.  If we are in a block of solids, move on.
        int *lProcRankForSite = mProcessorsForEachBlock[lBlockNumber].ProcessorRankForEachBlockSite;

        if (lProcRankForSite == NULL)
        {
          continue;
        }

        int m = -1;
        int partial_visited_fluid_sites = 0;

        for (int site_i = lBlockCoordI * block_size; site_i < lBlockCoordI
            * block_size + block_size; site_i++)
        {
          for (int site_j = lBlockCoordJ * block_size; site_j < lBlockCoordJ
              * block_size + block_size; site_j++)
          {
            for (int site_k = lBlockCoordK * block_size; site_k < lBlockCoordK
                * block_size + block_size; site_k++)
            {
              // Move on if the site is solid (ProcessorRankForEachBlockSite = 1 << 30) or has already been assigned
              // to a rank (0 <= ProcessorRankForEachBlockSite < 1 << 30).  ProcessorRankForEachBlockSite is allocated and initialised
              // in lbmReadConfig in io.cc.
              if (lProcRankForSite[++m] != marker)
              {
                continue;
              }
              // We have found an unvisited fluid site to start growing the subdomain from.
              // Assign it to the rank and update the fluid site counters.
              lProcRankForSite[m] = proc_count;

              if (IsCurrentProcRank(proc_count))
              {
                ++my_sites;
              }
              ++partial_visited_fluid_sites;

              if (unitLevel == 0)
              {
                ++mFluidSitesOnEachProcessor[proc_count];
              }

              // Sites on the edge of the mClusters at the start of the current graph growing partitioning step.
              int sites_a = 1;

              // Record the location of this initial site.
              SiteLocation *site_location_a_p = &site_location_a[0];
              site_location_a_p->i = site_i;
              site_location_a_p->j = site_j;
              site_location_a_p->k = site_k;

              // The subdomain can grow.
              bool are_fluid_sites_incrementing = true;

              // While the region can grow (i.e. it is not bounded by solids or visited
              // sites), and we need more sites on this particular rank.
              while (partial_visited_fluid_sites < fluid_sites_per_unit
                  && are_fluid_sites_incrementing)
              {
                // Sites added to the edge of the mClusters during the iteration.
                int sites_b = 0;
                are_fluid_sites_incrementing = false;

                // For sites on the edge of the domain (sites_a), deal with the neighbours.
                for (int index_a = 0; index_a < sites_a
                    && partial_visited_fluid_sites < fluid_sites_per_unit; index_a++)
                {
                  site_location_a_p = &site_location_a[index_a];

                  for (unsigned int l = 1; l < D3Q15::NUMVECTORS
                      && partial_visited_fluid_sites < fluid_sites_per_unit; l++)
                  {
                    // Record neighbour location.
                    int neigh_i = site_location_a_p->i + D3Q15::CX[l];
                    int neigh_j = site_location_a_p->j + D3Q15::CY[l];
                    int neigh_k = site_location_a_p->k + D3Q15::CZ[l];

                    // Move on if neighbour is outside the bounding box.
                    if (neigh_i == -1 || neigh_i == sites_x)
                      continue;
                    if (neigh_j == -1 || neigh_j == sites_y)
                      continue;
                    if (neigh_k == -1 || neigh_k == sites_z)
                      continue;

                    // Move on if the neighbour is in a block of solids (in which case
                    // the pointer to ProcessorRankForEachBlockSite is NULL) or it is solid or has already
                    // been assigned to a rank (in which case ProcessorRankForEachBlockSite != -1).  ProcessorRankForEachBlockSite
                    // was initialized in lbmReadConfig in io.cc.

                    // Pointer to the rank on which a particular fluid site
                    // resides.
                    int * proc_id_p = GetProcIdFromGlobalCoords(neigh_i,
                                                                neigh_j,
                                                                neigh_k);

                    if (proc_id_p == NULL || *proc_id_p != marker)
                    {
                      continue;
                    }
                    // Set the rank for a neighbour and update the fluid site counters.
                    *proc_id_p = proc_count;
                    ++partial_visited_fluid_sites;

                    if (unitLevel == 0)
                    {
                      ++mFluidSitesOnEachProcessor[proc_count];
                    }

                    if (IsCurrentProcRank(proc_count))
                    {
                      ++my_sites;
                    }
                    // Neighbour was found, so the region can grow.
                    are_fluid_sites_incrementing = true;

                    // If the new layer of neighbours is too large, allocate more
                    // memory.
                    if (sites_b == sites_buffer_size)
                    {
                      sites_buffer_size *= 2;
                      site_location_a
                          = (SiteLocation *) realloc(site_location_a,
                                                     sizeof(SiteLocation)
                                                         * sites_buffer_size);
                      site_location_b
                          = (SiteLocation *) realloc(site_location_b,
                                                     sizeof(SiteLocation)
                                                         * sites_buffer_size);
                    }

                    // Record the location of the neighbour.
                    SiteLocation * site_location_b_p =
                        &site_location_b[sites_b];
                    site_location_b_p->i = neigh_i;
                    site_location_b_p->j = neigh_j;
                    site_location_b_p->k = neigh_k;
                    ++sites_b;
                  }
                }
                // When the new layer of edge sites has been found, swap the buffers for
                // the current and new layers of edge sites.
                site_location_a_p = site_location_a;
                site_location_a = site_location_b;
                site_location_b = site_location_a_p;
                sites_a = sites_b;
              }
              // If we have enough sites, we have finished.
              if (partial_visited_fluid_sites >= fluid_sites_per_unit)
              {
                ++proc_count;
                if (unitLevel == 0)
                {
                  unvisited_fluid_sites -= partial_visited_fluid_sites;
                  fluid_sites_per_unit
                      = (int) ceil((double) unvisited_fluid_sites
                          / (double) (mProcessorCount - proc_count));
                }
                partial_visited_fluid_sites = 0;
              }
              // If not, we have to start growing a different region for the same rank:
              // region expansions could get trapped.
            }
          }
        }
      }
    }
  }

  delete[] site_location_b;
  delete[] site_location_a;
}

/*!
 This is called from the main function.  First function to deal with processors.
 The domain partitioning technique and the management of the
 buffers useful for the inter-processor communications are
 implemented in this function.  The domain decomposition is based
 on a graph growing partitioning technique.
 */
void Net::Initialise(int iTotalFluidSites)
{
  // Allocations.  fluid sites will store actual number of fluid
  // sites per proc.  Site location will store up to 10000 of some
  // sort of coordinate.
  mFluidSitesOnEachProcessor = new int[mProcessorCount];

  double seconds;
  int site_i, site_j, site_k; // Global coordinates of a site.
  int neigh_i, neigh_j, neigh_k; // Global coordinates of a neighbour site.
  int i, j, k; // Global block index.
  unsigned int l; // Index for neighbours of a site.
  int m; // Site index on a paricular block.
  int n; // Global block index.
  int mm; // Index of processors surrounding this one.
  int unvisited_fluid_sites; // Fluid sites not yet visited.
  int proc_count; // Rank we are looking at.
  int fluid_sites_per_unit; // Fluid sites per rank.
  int neigh_proc_index;
  int my_sites_p; // Sites residing on this rank.

  int is_inter_site, is_inner_site; // Whether the site is on the edge of the subdomain of a rank.
  int collision_offset[2][COLLISION_TYPES];
  int flag; // Whether a neighbouring process has been listed (0) or
  // not (1).

  short int *f_data_p;
  unsigned int *site_data; // Local variable to store data for sites on this rank.
  unsigned int *site_data_p;
  unsigned int site_map;
  bool *is_my_block; // Array to store whether any sites on a block are fluid
  // sites residing on this rank.

  DataBlock *data_block_p; // Pointer fitting into DataBlock type struct (a single-member struct).
  DataBlock *map_block_p; // Pointer fitting into DataBlock type struct (a single-member struct).
  ProcBlock *proc_block_p; // Pointer fitting into ProcBlock type struct (a single-member struct).
  NeighProc *neigh_proc_p; // Pointer fitting into NeighProc type struct.
  net_site_nor = NULL;
  net_site_data = NULL;
  cut_distances = NULL;

  my_sites = 0;

  // a fast graph growing partitioning technique which spans the data
  // set only once is implemented here; the data set is explored by
  // means of the arrays "site_location_a[]" and "site_location_b[]"

  // Find the maximum number of fluid sites per process.  If steering,
  // leave one process out.
#ifndef NO_STEER 
  if (mProcessorCount == 1)
  {
    fluid_sites_per_unit = (int) ceil((double) iTotalFluidSites
        / (double) mProcessorCount);
    proc_count = 0;
  }
  else
  {
    fluid_sites_per_unit = (int) ceil((double) iTotalFluidSites
        / (double) (mProcessorCount - 1));
    proc_count = 1;
  }
#else
  fluid_sites_per_unit = (int)ceil((double)iTotalFluidSites / (double)mProcessorCount);
  proc_count = 0;
#endif

  for (n = 0; n < mProcessorCount; n++)
  {
    mFluidSitesOnEachProcessor[n] = 0;
  }

  unvisited_fluid_sites = iTotalFluidSites;

  seconds = hemelb::util::myClock();

  if (net_machines == 1 || net_machines == mProcessorCount) // If one machine or one machine per proc.
  {
    ThisNeedsRenaming(proc_count, fluid_sites_per_unit, unvisited_fluid_sites,
                      -1, 0);
  }
  else
  {
    for (int unit_level = 1; unit_level >= 0; unit_level--)
    {
      int up_units_max;

      if (unit_level == 1)
      {
        up_units_max = 1;
        fluid_sites_per_unit = 0;
      }
      else
      {
        up_units_max = net_machines;

        fluid_sites_per_unit = (int) ceil((double) unvisited_fluid_sites
            / (double) (mProcessorCount - 1));
        proc_count = 1;

      }
      for (int up_unit = 0; up_unit < up_units_max; up_unit++)
      {
        int marker;

        if (unit_level == 1)
        {
          marker = -1;

          proc_count = mProcessorCount;

          int machine_id = proc_count - mProcessorCount;

          double weight;

          weight = (double) (procs_per_machine[machine_id] * mProcessorCount)
              / (double) (mProcessorCount - 1);

          fluid_sites_per_unit = (int) ceil((double) iTotalFluidSites * weight
              / net_machines);
        }
        else
        {
          marker = mProcessorCount + up_unit;
        }

        ThisNeedsRenaming(proc_count, fluid_sites_per_unit,
                          unvisited_fluid_sites, marker, unit_level);
      }
    }
  }

  dd_time = hemelb::util::myClock() - seconds;
  seconds = hemelb::util::myClock();

  // A map between the two-level data representation and the 1D
  // compact one is created here.

  // Allocate an array of structures to store the fluid site identifiers for each block.
  map_block = new DataBlock[blocks];

  for (n = 0; n < blocks; n++)
  {
    map_block[n].site_data = NULL;
  }

  // Local array to store site data for this rank.
  site_data = new unsigned int[my_sites];

  // Allocate blocks.
  is_my_block = new bool[blocks];

  for (n = 0; n < blocks; n++)
  {
    is_my_block[n] = 0;
  }
  my_sites_p = 0;

  for (n = 0; n < blocks; n++)
  {
    // If we are in a block of solids, move to the next block.
    data_block_p = &data_block[n];

    if (data_block_p->site_data == NULL)
    {
      continue;
    }
    // If we have some fluid sites, point to mProcessorsForEachBlock and map_block.
    proc_block_p = &mProcessorsForEachBlock[n];
    map_block_p = &map_block[n];
    map_block_p->site_data = new unsigned int[sites_in_a_block];

    // map_block[n].site_data is set to the fluid site identifier on this rank or (1U << 31U) if a site is solid
    // or not on this rank.  site_data is indexed by fluid site identifier and set to the site_data.
    for (m = 0; m < sites_in_a_block; m++)
    {
      if (IsCurrentProcRank(proc_block_p->ProcessorRankForEachBlockSite[m]))
      {
        if ( (data_block_p->site_data[m] & SITE_TYPE_MASK) != SOLID_TYPE)
        {
          map_block_p->site_data[m] = my_sites_p;
          site_data[my_sites_p] = data_block_p->site_data[m];
          ++my_sites_p;
        }
        else
        {
          map_block_p->site_data[m] = (1U << 31U);
        }
        is_my_block[n] = 1;
      }
      else
      {
        map_block_p->site_data[m] = (1U << 31U);
      }
    }
  }

  // Free data_block.
  for (n = 0; n < blocks; n++)
  {
    if (data_block[n].site_data != NULL)
    {
      delete[] data_block[n].site_data;
      data_block[n].site_data = NULL;
    }
  }
  delete[] data_block;
  data_block = NULL;

  // If we are in a block of solids, we set map_block[n].site_data to NULL.
  for (n = 0; n < blocks; n++)
  {
    if (is_my_block[n])
      continue;

    delete[] map_block[n].site_data;
    map_block[n].site_data = NULL;

    if (wall_block[n].wall_data != NULL)
    {
      delete[] wall_block[n].wall_data;
      wall_block[n].wall_data = NULL;
    }
  }
  delete[] is_my_block;

  // The numbers of inter- and intra-machine neighbouring processors,
  // interface-dependent and independent fluid sites and shared
  // distribution functions of the reference processor are calculated
  // here.  neigh_proc is a static array that is declared in config.h.

  neigh_procs = 0;

  for (m = 0; m < NEIGHBOUR_PROCS_MAX; m++)
  {
    neigh_proc[m].fs = 0; // fs within NeighProc struct within the net struct.
  }

  my_inter_sites = 0;
  my_inner_sites = 0;

  for (m = 0; m < COLLISION_TYPES; m++)
  {
    my_inter_collisions[m] = 0;
    my_inner_collisions[m] = 0;
  }
  shared_fs = 0; // shared fs within Net struct.

  my_sites_p = 0;

  n = -1;

  // Here, i, j and k are not the block coordinates, but the block coords * block_size.
  for (i = 0; i < sites_x; i += block_size)
    for (j = 0; j < sites_y; j += block_size)
      for (k = 0; k < sites_z; k += block_size)
      {
        map_block_p = &map_block[++n];

        if (map_block_p->site_data == NULL)
          continue;

        proc_block_p = &mProcessorsForEachBlock[n];

        m = -1;

        for (site_i = i; site_i < i + block_size; site_i++)
          for (site_j = j; site_j < j + block_size; site_j++)
            for (site_k = k; site_k < k + block_size; site_k++)
            {
              m++;
              if (!IsCurrentProcRank(
                                     proc_block_p->ProcessorRankForEachBlockSite[m]))
              {
                continue;
              }
              is_inter_site = 0;
              is_inner_site = 1;

              for (l = 1; l < D3Q15::NUMVECTORS; l++)
              {
                neigh_i = site_i + D3Q15::CX[l];
                neigh_j = site_j + D3Q15::CY[l];
                neigh_k = site_k + D3Q15::CZ[l];

                int *proc_id_p = GetProcIdFromGlobalCoords(neigh_i, neigh_j,
                                                           neigh_k);

                // Move on if the neighbour is in a block of solids (in which case
                // the pointer to ProcessorRankForEachBlockSite is NULL) or it is solid (in which case ProcessorRankForEachBlockSite ==
                // 1 << 30) or the neighbour is also on this rank.  ProcessorRankForEachBlockSite was initialized
                // in lbmReadConfig in io.cc.
                if (proc_id_p == NULL || IsCurrentProcRank(*proc_id_p)
                    || *proc_id_p == (1 << 30))
                {
                  continue;
                }
                is_inner_site = 0;
                is_inter_site = 1;

                // The first time, we set mm = 0, flag = 1, but net_neigh_procs = 0, so
                // the loop is not executed.
                for (mm = 0, flag = 1; mm < neigh_procs && flag; mm++)
                {
                  // Check whether the rank for a particular neighbour has already been
                  // used for this processor.  If it has, set flag to zero.
                  neigh_proc_p = &neigh_proc[mm];

                  // If ProcessorRankForEachBlockSite is equal to a neigh_proc that has alredy been listed.
                  if (*proc_id_p == neigh_proc_p->id)
                  {
                    flag = 0;
                    ++neigh_proc_p->fs;
                    ++shared_fs;
                  }
                }
                // If flag is 1, we need a new neighbouring processor.
                if (flag)
                {
                  if (neigh_procs == NEIGHBOUR_PROCS_MAX)
                  {
                    printf(
                           " too many intra machine, inter processor neighbours\n");
                    printf(" the execution is terminated\n");
                    Abort();
                  }
                  // Store rank of neighbour in >neigh_proc[neigh_procs]
                  neigh_proc_p = &neigh_proc[neigh_procs];
                  neigh_proc_p->id = *proc_id_p;
                  ++neigh_proc_p->fs;
                  ++neigh_procs;
                  ++shared_fs;
                }
              }
              // Collision Type set here. map_block site data is renumbered according to
              // fluid site numbers within a particular collision type.
              if (GetCollisionType(site_data[my_sites_p]) == FLUID)
              {
                l = 0;
              }
              else if (GetCollisionType(site_data[my_sites_p]) == EDGE)
              {
                l = 1;
              }
              else if (GetCollisionType(site_data[my_sites_p]) == INLET)
              {
                l = 2;
              }
              else if (GetCollisionType(site_data[my_sites_p]) == OUTLET)
              {
                l = 3;
              }
              else if (GetCollisionType(site_data[my_sites_p])
                  == (INLET | EDGE))
              {
                l = 4;
              }
              else if (GetCollisionType(site_data[my_sites_p]) == (OUTLET
                  | EDGE))
              {
                l = 5;
              }

              ++my_sites_p;

              if (is_inner_site)
              {
                ++my_inner_sites;

                if (l == 0)
                {
                  map_block_p->site_data[m] = my_inner_collisions[l];
                }
                else
                {
                  map_block_p->site_data[m] = 50000000 * (10 + (l - 1))
                      + my_inner_collisions[l];
                }
                ++my_inner_collisions[l];
              }
              else if (is_inter_site)
              {
                ++my_inter_sites;

                if (l == 0)
                {
                  map_block_p->site_data[m] = 1000000000
                      + my_inter_collisions[l];
                }
                else
                {
                  map_block_p->site_data[m] = 50000000 * (20 + l)
                      + my_inter_collisions[l];
                }
                ++my_inter_collisions[l];
              }
            }
      }
  // Calculte the number of each type of collision.
  collision_offset[0][0] = 0;

  for (l = 1; l < COLLISION_TYPES; l++)
  {
    collision_offset[0][l] = collision_offset[0][l - 1] + my_inner_collisions[l
        - 1];
  }
  collision_offset[1][0] = my_inner_sites;

  for (l = 1; l < COLLISION_TYPES; l++)
  {
    collision_offset[1][l] = collision_offset[1][l - 1] + my_inter_collisions[l
        - 1];
  }
  for (n = 0; n < blocks; n++)
  {
    map_block_p = &map_block[n];

    // If we are in a block of solids, continue.
    if (map_block_p->site_data == NULL)
      continue;

    for (m = 0; m < sites_in_a_block; m++)
    {
      site_data_p = &map_block_p->site_data[m];

      // If the site is solid, continue.
      if (*site_data_p & (1U << 31U))
        continue;

      // 0th collision type for inner sites, so don't do anything.
      if (*site_data_p < 500000000)
      {
        continue;
      }
      // Renumber the sites in map_block so that the numbers are compacted together.  We have
      // collision offset to tell us when one collision type ends and another starts.
      for (l = 1; l < COLLISION_TYPES; l++)
      {
        if (*site_data_p >= 50000000 * (10 + (l - 1)) && *site_data_p
            < 50000000 * (10 + l))
        {
          *site_data_p += collision_offset[0][l] - 50000000 * (10 + (l - 1));
          break;
        }
      }
      for (l = 0; l < COLLISION_TYPES; l++)
      {
        if (*site_data_p >= 50000000 * (20 + l) && *site_data_p < 50000000
            * (20 + (l + 1)))
        {
          *site_data_p += collision_offset[1][l] - 50000000 * (20 + l);
          break;
        }
      }
    }
  }

  // Allocate f_old and f_new according to the number of sites on the process.  The extra site
  // is there for when we would stream into a solid site during the simulation, which avoids
  // an if condition at every timestep at every boundary site.  We also allocate space for the
  // shared distribution functions.  We need twice as much space when we check the convergence
  // and the extra distribution functions are
  f_old = new double[my_sites * D3Q15::NUMVECTORS + 1 + shared_fs];
  f_new = new double[my_sites * D3Q15::NUMVECTORS + 1 + shared_fs];

  // the precise interface-dependent data (interface-dependent fluid
  // site locations and identifiers of the distribution functions
  // streamed between different partitions) are collected and the
  // buffers needed for the communications are set from here

  f_data = new short int[4 * shared_fs];

  // Allocate the index in which to put the distribution functions received from the other
  // process.
  f_recv_iv = new int[shared_fs];

  // Reset to zero again.
  shared_fs = 0;

  for (n = 0; n < neigh_procs; n++)
  {
    // f_data compacted according to number of shared f_s on each process.
    // f_data will be set later.
    neigh_proc[n].f_data = &f_data[shared_fs << 2];

    // Pointing to a few things, but not setting any variables.
    // f_head points to start of shared_fs.
    neigh_proc[n].f_head = my_sites * D3Q15::NUMVECTORS + 1 + shared_fs;

    neigh_proc[n].f_recv_iv = &f_recv_iv[shared_fs];

    shared_fs += neigh_proc[n].fs;
    neigh_proc[n].fs = 0;// This is set back to 0.
  }
  if (my_sites > 0)
  {
    // f_id is allocated so we know which sites to get information from.
    f_id = new int[my_sites * D3Q15::NUMVECTORS];

    net_site_data = new unsigned int[my_sites];

    if (lbm_stress_type == SHEAR_STRESS)
    {
      net_site_nor = new double[my_sites * 3];
      cut_distances = new double[my_sites * (D3Q15::NUMVECTORS - 1)];
    }
  }
  from_proc_id_to_neigh_proc_index = new short int[mProcessorCount];

  for (m = 0; m < mProcessorCount; m++)
  {
    from_proc_id_to_neigh_proc_index[m] = -1;
  }
  // Get neigh_proc_index from ProcessorRankForEachBlockSite.
  for (m = 0; m < neigh_procs; m++)
  {
    from_proc_id_to_neigh_proc_index[neigh_proc[m].id] = m;
  }
  my_sites_p = 0;

  n = -1;

  for (i = 0; i < sites_x; i += block_size)
    for (j = 0; j < sites_y; j += block_size)
      for (k = 0; k < sites_z; k += block_size)
      {
        map_block_p = &map_block[++n];

        if (map_block_p->site_data == NULL)
          continue;

        proc_block_p = &mProcessorsForEachBlock[n];

        m = -1;

        for (site_i = i; site_i < i + block_size; site_i++)
          for (site_j = j; site_j < j + block_size; site_j++)
            for (site_k = k; site_k < k + block_size; site_k++)
            {
              // If a site is not on this process, continue.
              m++;
              if (!IsCurrentProcRank(
                                     proc_block_p->ProcessorRankForEachBlockSite[m]))
              {
                continue;
              }
              // Get site data, which is the number of the fluid site on this proc..
              site_map = map_block_p->site_data[m];

              // set f_id.
              f_id[site_map * D3Q15::NUMVECTORS + 0] = site_map
                  * D3Q15::NUMVECTORS + 0;

              for (l = 1; l < D3Q15::NUMVECTORS; l++)
              {
                // Work out positions of neighbours.
                neigh_i = site_i + D3Q15::CX[l];
                neigh_j = site_j + D3Q15::CY[l];
                neigh_k = site_k + D3Q15::CZ[l];

                // initialize f_id to the rubbish site.
                f_id[site_map * D3Q15::NUMVECTORS + l] = my_sites
                    * D3Q15::NUMVECTORS;

                // You know which process the neighbour is on.
                int *proc_id_p = GetProcIdFromGlobalCoords(neigh_i, neigh_j,
                                                           neigh_k);

                if (proc_id_p == NULL || *proc_id_p == 1 << 30)
                {
                  continue;
                }
                // Pointer to the neihgbour.
                site_data_p = netSiteMapPointer(neigh_i, neigh_j, neigh_k);

                // If on the same proc, set f_id of the
                // current site and direction to the
                // site and direction that it sends to.
                // If we check convergence, the data for
                // each site is split into that for the
                // current and previous cycles.
                if (IsCurrentProcRank(*proc_id_p))
                {
                  f_id[site_map * D3Q15::NUMVECTORS + l] = *site_data_p
                      * D3Q15::NUMVECTORS + l;

                  continue;
                }
                neigh_proc_index = from_proc_id_to_neigh_proc_index[*proc_id_p];

                // You have neigh proc again.
                neigh_proc_p = &neigh_proc[neigh_proc_index];

                // This stores some coordinates.  We
                // still need to know the site number.
                // neigh_proc[ n ].f_data is now
                // set as well, since this points to
                // f_data.  Every process has data for
                // its neighbours which say which sites
                // on this process are shared with the
                // neighbour.
                f_data_p = &neigh_proc_p->f_data[neigh_proc_p->fs << 2];
                f_data_p[0] = site_i;
                f_data_p[1] = site_j;
                f_data_p[2] = site_k;
                f_data_p[3] = l;
                ++neigh_proc_p->fs; // We recount this again.
              }
              // This is used in Calculate BC in IO.
              net_site_data[site_map] = site_data[my_sites_p];

              if (lbm_stress_type == SHEAR_STRESS)
              {
                if (GetCollisionType(net_site_data[site_map]) & EDGE)
                {
                  for (l = 0; l < 3; l++)
                    net_site_nor[site_map * 3 + l]
                        = wall_block[n].wall_data[m].wall_nor[l];

                  for (l = 0; l < (D3Q15::NUMVECTORS - 1); l++)
                    cut_distances[site_map * (D3Q15::NUMVECTORS - 1) + l]
                        = wall_block[n].wall_data[m].cut_dist[l];
                }
                else
                {
                  net_site_nor[site_map * 3] = 1.0e+30;
                }
              }
              ++my_sites_p;
            }
      }
  delete[] site_data;

  // point-to-point communications are performed to match data to be
  // sent to/receive from different partitions; in this way, the
  // communication of the locations of the interface-dependent fluid
  // sites and the identifiers of the distribution functions which
  // propagate to different partitions is avoided (only their values
  // will be communicated). It's here!

  // Allocate the request variable.
#ifndef NOMPI
  req = new MPI_Request*[COMMS_LEVELS];

  for (m = 0; m < COMMS_LEVELS; m++)
  {
    req[m] = new MPI_Request[2 * mProcessorCount];
  }
#endif

  for (m = 0; m < neigh_procs; m++)
  {
    neigh_proc_p = &neigh_proc[m];

    // One way send receive.  The lower numbered mProcessorCount send and the higher numbered ones receive.
    // It seems that, for each pair of processors, the lower numbered one ends up with its own
    // edge sites and directions stored and the higher numbered one ends up with those on the
    // other processor.
#ifndef NOMPI
    if (neigh_proc_p->id > mRank)
    {
      err = MPI_Isend(&neigh_proc_p->f_data[0], neigh_proc_p->fs * 4,
                      MPI_SHORT, neigh_proc_p->id, 10, MPI_COMM_WORLD,
                      &req[0][m]);
    }
    else
    {
      err = MPI_Irecv(&neigh_proc_p->f_data[0], neigh_proc_p->fs * 4,
                      MPI_SHORT, neigh_proc_p->id, 10, MPI_COMM_WORLD,
                      &req[0][neigh_procs + m]);
    }
#endif
  }
  for (m = 0; m < neigh_procs; m++)
  {
    neigh_proc_p = &neigh_proc[m];

    if (neigh_proc_p->id > mRank)
    {
#ifndef NOMPI
      err = MPI_Wait(&req[0][m], status);
#endif
    }
    else
    {
#ifndef NOMPI
      err = MPI_Wait(&req[0][neigh_procs + m], status);
#endif
      // Now we sort the situation so that each process has its own sites.
      for (n = 0; n < neigh_proc_p->fs * 4; n += 4)
      {
        f_data_p = &neigh_proc_p->f_data[n];

        l = f_data_p[3];
        f_data_p[0] += D3Q15::CX[l];
        f_data_p[1] += D3Q15::CY[l];
        f_data_p[2] += D3Q15::CZ[l];
        f_data_p[3] = D3Q15::INVERSEDIRECTIONS[l];
      }
    }
  }

  int f_count = my_sites * D3Q15::NUMVECTORS;

  for (m = 0; m < neigh_procs; m++)
  {
    neigh_proc_p = &neigh_proc[m];

    for (n = 0; n < neigh_proc_p->fs; n++)
    {
      // Get coordinates and direction of the distribution function to be sent to another process.
      f_data_p = &neigh_proc_p->f_data[n * 4];
      i = f_data_p[0];
      j = f_data_p[1];
      k = f_data_p[2];
      l = f_data_p[3];

      // Get the fluid site number of site that will send data to another process.
      site_map = *netSiteMapPointer(i, j, k);

      // Set f_id to the element in the send buffer that we put the updated
      // distribution functions in.
      f_id[site_map * D3Q15::NUMVECTORS + l] = ++f_count;

      // Set the place where we put the received distribution functions, which is
      // f_new[number of fluid site that sends, inverse direction].
      neigh_proc_p->f_recv_iv[n] = site_map * D3Q15::NUMVECTORS
          + D3Q15::INVERSEDIRECTIONS[l];
    }
  }
  // neigh_prc->f_data was only set as a pointer to f_data, not allocated.  In this line, we 
  // are freeing both of those.
  delete[] f_data;

  bm_time = hemelb::util::myClock() - seconds;
}

void Net::ReceiveFromNeighbouringProcessors()
{
#ifndef NOMPI
  for (int m = 0; m < neigh_procs; m++)
  {
    NeighProc *neigh_proc_p = &neigh_proc[m];

    err = MPI_Irecv(&f_old[neigh_proc_p->f_head], neigh_proc_p->fs, MPI_DOUBLE,
                    neigh_proc_p->id, 10, MPI_COMM_WORLD, &req[0][m]);
  }
#endif
}

int Net::GetMachineCount()
{
  return net_machines;
}

void Net::SendToNeighbouringProcessors()
{
#ifndef NOMPI
  for (int m = 0; m < neigh_procs; m++)
  {
    NeighProc *neigh_proc_p = &neigh_proc[m];

    err = MPI_Isend(&f_new[neigh_proc_p->f_head], neigh_proc_p->fs, MPI_DOUBLE,
                    neigh_proc_p->id, 10, MPI_COMM_WORLD, &req[0][neigh_procs
                        + m]);
  }
#endif
}

void Net::UseDataFromNeighbouringProcs()
{
#ifndef NOMPI
  for (int m = 0; m < neigh_procs; m++)
  {
    err = MPI_Wait(&req[0][m], status);
    err = MPI_Wait(&req[0][neigh_procs + m], status);
  }
#endif

  // Copy the distribution functions received from the neighbouring
  // processors into the destination buffer "f_new".
  for (int i = 0; i < shared_fs; i++)
  {
    f_new[f_recv_iv[i]] = f_old[neigh_proc[0].f_head + i];
  }
}

Net::Net(int &iArgumentCount, char* iArgumentList[])
{
#ifndef NOMPI
  int thread_level_provided;

  err = MPI_Init_thread(&iArgumentCount, &iArgumentList, MPI_THREAD_FUNNELED,
                        &thread_level_provided);
  err = MPI_Comm_size(MPI_COMM_WORLD, &mProcessorCount);
  err = MPI_Comm_rank(MPI_COMM_WORLD, &mRank);

  if (IsCurrentProcTheIOProc())
  {
    printf("thread_level_provided %i\n", thread_level_provided);
  }
#else
  mProcessorCount = 1;
  mRank = 0;
#endif
}

/*!
 Free the allocated data.
 */
Net::~Net()
{

#ifndef NOMPI
  err = MPI_Finalize();
#endif

  delete[] from_proc_id_to_neigh_proc_index;

  delete[] f_recv_iv;

  if (lbm_stress_type == SHEAR_STRESS && my_sites > 0)
  {
    delete[] net_site_nor;
    delete[] cut_distances;

    for (int i = 0; i < blocks; i++)
    {
      if (wall_block[i].wall_data != NULL)
      {
        delete[] wall_block[i].wall_data;
      }
    }
    delete[] wall_block;
  }

  for (int i = 0; i < blocks; i++)
  {
    if (map_block[i].site_data != NULL)
    {
      delete[] map_block[i].site_data;
    }
  }
  delete[] map_block;

  for (int i = 0; i < blocks; i++)
  {
    if (mProcessorsForEachBlock[i].ProcessorRankForEachBlockSite != NULL)
    {
      delete[] mProcessorsForEachBlock[i].ProcessorRankForEachBlockSite;
    }
  }
  delete[] mProcessorsForEachBlock;

  if (my_sites > 0)
  {
    delete[] net_site_data;
    delete[] f_id;
  }
  delete[] f_new;
  delete[] f_old;

#ifndef NOMPI
  for (int i = 0; i < COMMS_LEVELS; i++)
    delete[] req[i];
  delete[] req;
#endif

  delete[] mFluidSitesOnEachProcessor;

  delete[] procs_per_machine;
  delete[] machine_id;
}

