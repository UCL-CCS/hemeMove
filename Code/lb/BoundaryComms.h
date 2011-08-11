#ifndef HEMELB_LB_BOUNDARYCOMMS_H
#define HEMELB_LB_BOUNDARYCOMMS_H

#include "net/IteratedAction.h"
#include "geometry/LatticeData.h"
#include "lb/SimulationState.h"
#include "util/UnitConverter.h"

namespace hemelb
{
  namespace lb
  {

    /*
     * This class deals with in/outlet boundary conditions. This class is not concerned with updating
     * the densities. It assumes that it has already done when Broadcast is called. Every in/outlet will
     * have a corresponding communicator which includes the BCproc and all processes that contain the
     * given in/outlet. Once per cycle BroadcastBoundaryDensities needs to be called so that BCproc
     * updates the relevant processes with the new values.
     */
    class BoundaryComms : public net::IteratedAction
    {
      public:
        BoundaryComms(const geometry::LatticeData* iLatDat,
                      SimConfig* iSimConfig,
                      SimulationState* iSimState,
                      util::UnitConverter* iUnits);
        ~BoundaryComms();

        void RequestComms();
        void EndIteration();
        void Reset();

        void WaitForComms(const int index, unsigned int IOtype);

        void InitialiseBoundaryDensities();
        void SendBoundaryDensities();
        void CalculateBC(distribn_t f[],
                         hemelb::geometry::LatticeData::SiteType iSiteType,
                         unsigned int iBoundaryId,
                         distribn_t *density,
                         distribn_t *vx,
                         distribn_t *vy,
                         distribn_t *vz,
                         distribn_t f_neq[]);

        distribn_t GetInletDensity(int i);
        distribn_t GetOutletDensity(int i);

        // The densities
        // Unfortunately made public, because of a lot of old code relying on this being in LBM
        // TODO: Make private
        distribn_t *inlet_density, *outlet_density;
        distribn_t *inlet_density_avg, *outlet_density_avg;
        distribn_t *inlet_density_amp, *outlet_density_amp;

      private:
        proc_t BCproc; // Process responsible for sending out BC info
        distribn_t *inlet_density_cycle, *outlet_density_cycle;

        void ReadParameters();
        void allocateInlets();
        void allocateOutlets();

        // Total number of inlets/outlets in simulation
        int nTotInlets;
        int nTotOutlets;

        distribn_t *inlet_density_phs, *outlet_density_phs;

        // Number of inlets/outlets on this process
        int nInlets;
        int nOutlets;

        // List of indices of inlets/outlets on this process
        std::vector<int> inlets;
        std::vector<int> outlets;

        // These are only assigned on the BC proc as it is the only one that needs to know
        // which proc has which IOlet
        int* inletRequestOffset;
        int* outletRequestOffset;
        int* nInletProcs;
        int* nOutletProcs;
        int** inletProcsList;
        int** outletProcsList;

        MPI_Request* inlet_request;
        MPI_Request* outlet_request;
        MPI_Status* inlet_status;
        MPI_Status* outlet_status;

        bool* inletCommAlreadyChecked;
        bool* outletCommAlreadyChecked;

        SimulationState* mState;
        SimConfig* mSimConfig;
        util::UnitConverter* mUnits;
    };

  }
}

#endif /* HEMELB_LB_BOUNDARYCOMMS_H */
