#ifndef HEMELB_NET_H
#define HEMELB_NET_H

#include "constants.h"
#include "mpiInclude.h"
#include "D3Q15.h"
#include "SimConfig.h"

#include "lb/GlobalLatticeData.h"
#include "lb/LocalLatticeData.h"
// Superficial site data
struct WallData
{
    // estimated boundary normal (if the site is an inlet/outlet site)
    double boundary_nor[3];
    // estimated minimum distance (in lattice units) from the
    // inlet/outlet boundary;
    double boundary_dist;
    // estimated wall normal (if the site is close to the wall);
    double wall_nor[3];
    // estimated minimum distance (in lattice units) from the wall;
    // if the site is close to the wall surface
    double wall_dist;
    // cut distances along the 14 non-zero lattice vectors;
    // each one is between 0 and 1 if the surface cuts the corresponding
    // vector or is equal to 1e+30 otherWallDatawise
    double cut_dist[D3Q15::NUMVECTORS - 1];

};

// For each global block,
// *ProcessorRankForEachBlockSite is an array containing the ranks on which individual
// lattice sites reside.
// WallBlock is a member of the structure Net and is employed to store the data
// regarding the wall, inlet and outlet sites.
// Block means macrocell of fluid sites (voxels).
// site_data[] is an array containing individual lattice site data
// within a global block.
struct DataBlock
{
    int *ProcessorRankForEachBlockSite;
    WallData *wall_data;
    unsigned int *site_data;
};

class Net
{
  public:
    Net(int &iArgumentCount, char* iArgumentList[]);
    ~Net();

    void Abort();
    bool IsCurrentProcRank(int iRank);
    bool IsCurrentProcTheIOProc();
    char *GetCurrentProcIdentifier();

    // declarations of all the functions used
    int *GetProcIdFromGlobalCoords(int iSiteI,
                                   int iSiteJ,
                                   int iSiteK,
                                   hemelb::lb::GlobalLatticeData &iGlobLatDat);
    int netFindTopology();
    void Initialise(int iTotalFluidSites,
                    hemelb::lb::GlobalLatticeData &iGlobLatDat,
                    hemelb::lb::LocalLatticeData &iLocalLatDat);

    double GetCutDistance(int iSiteIndex, int iDirection) const;
    bool HasBoundary(int iSiteIndex, int iDirection) const;
    int GetBoundaryId(int iSiteIndex) const;
    double* GetNormalToWall(int iSiteIndex) const;

    void ReceiveFromNeighbouringProcessors(hemelb::lb::LocalLatticeData &bLocalLatDat);
    void SendToNeighbouringProcessors(hemelb::lb::LocalLatticeData &bLocalLatDat);
    void UseDataFromNeighbouringProcs(hemelb::lb::LocalLatticeData &bLocalLatDat);

    int GetMachineCount();

    // Currently needed for the destructor. TODO should probably do via vector or similar.
    int block_count;

    int mProcessorCount; // Number of processors.
    int err;
    int my_inner_sites; // Site on this process that do and do not need
    // information from neighbouring processors.
    int my_inner_collisions[COLLISION_TYPES]; // Number of collisions that only use data on this rank.
    int my_inter_collisions[COLLISION_TYPES]; // Number of collisions that require information from
    // other processors.
    int my_sites; // Number of fluid sites on this rank.
    int *mFluidSitesOnEachProcessor; // Array containing numbers of fluid sites on
    // each process.

    DataBlock *map_block; // See comment next to struct DataBlock.

    unsigned int GetCollisionType(unsigned int site_data);
    MPI_Status status[4];
    double dd_time, bm_time, fr_time;
    unsigned int *net_site_data;

    int GetDepths();
  private:
    // NeighProc is part of the Net (defined later in this file).  This object is an element of an array
    // (called neigh_proc[]) and comprises information about the neighbouring processes to this process.
    struct NeighProc
    {
        int id; // Rank of the neighbouring processor.
        int fs; // Number of distributions shared with neighbouring
        // processors.
        short int *f_data; // Coordinates of a fluid site that streams to the on
        // neighbouring processor "id" and
        // streaming direction
        int f_head;
        int *f_recv_iv;
    };
    // Some sort of coordinates.
    struct SiteLocation
    {
        short int i, j, k;
    };
    int *f_recv_iv;
    int my_inter_sites;
    unsigned int *netSiteMapPointer(int site_i,
                                    int site_j,
                                    int site_k,
                                    hemelb::lb::GlobalLatticeData &iGlobLatDat);
    void
    AssignFluidSitesToProcessors(int & proc_count,
                                 int & fluid_sites_per_unit,
                                 int & unvisited_fluid_sites,
                                 const int marker,
                                 const int unitLevel,
                                 hemelb::lb::GlobalLatticeData &iGlobLatDat);
    NeighProc neigh_proc[NEIGHBOUR_PROCS_MAX]; // See comment next to struct NeighProc.
    int shared_fs; // Number of distributions shared with neighbouring
    // processors.
    double *cut_distances;
    hemelb::SimConfig* mSimConfig;
    double *net_site_nor;
    int *machine_id;
    int *procs_per_machine;
    short int *from_proc_id_to_neigh_proc_index; // Turns ProcessorRankForEachBlockSite to neigh_proc_iindex.
    int neigh_procs; // Number of neighbouring processors.
    int mRank;
    MPI_Request **req;
    // 3 buffers needed for convergence-enabled simulations
    short int *f_data;
    int depths;
    int mMachineCount;
};

#endif // HEMELB_NET_H
