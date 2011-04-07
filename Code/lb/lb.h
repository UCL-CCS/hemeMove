#ifndef HEMELB_LB_H
#define HEMELB_LB_H

#include "net/net.h"
#include "net/IteratedAction.h"
#include "topology/NetworkTopology.h"
#include "lb/SimulationState.h"
#include "lb/collisions/Collisions.h"
#include "vis/ColPixel.h"
#include "SimConfig.h"

namespace hemelb
{
  namespace lb
  {
    class LBM : public net::IteratedAction
    {
      public:
        int total_fluid_sites;
        int inlets, outlets;
        unsigned int period;

        double ConvertPressureToLatticeUnits(double pressure) const;
        double ConvertPressureToPhysicalUnits(double pressure) const;
        double ConvertVelocityToLatticeUnits(double velocity) const;
        double ConvertStressToLatticeUnits(double stress) const;
        double ConvertStressToPhysicalUnits(double stress) const;
        double ConvertVelocityToPhysicalUnits(double velocity) const;

        LBM(hemelb::SimConfig *iSimulationConfig,
            net::Net* net,
            geometry::LatticeData* latDat,
            SimulationState* simState,
            const hemelb::topology::NetworkTopology * iNetTop);
        void Reset();
        ~LBM();

        void CalculateFlowFieldValues();
        void RecalculateTauViscosityOmega();
        void UpdateBoundaryDensities(int cycle_id, int time_step);
        void UpdateInletVelocities(int time_step);

        void Initialise(int* iFTranslator, vis::Control* iControl);

        void SetInitialConditions();

        void
        WriteConfigParallel(hemelb::lb::Stability stability, std::string output_file_name);

        double GetMinPhysicalPressure();
        double GetMaxPhysicalPressure();
        double GetMinPhysicalVelocity();
        double GetMaxPhysicalVelocity();
        double GetMinPhysicalStress();
        double GetMaxPhysicalStress();

        void InitMinMaxValues(void);

        double GetAverageInletVelocity(int iInletNumber);
        double GetPeakInletVelocity(int iInletNumber);

        void ReadVisParameters();

        void CalculateMouseFlowField(float densityIn,
                                     float stressIn,
                                     double &mouse_pressure,
                                     double &mouse_stress,
                                     double density_threshold_min,
                                     double density_threshold_minmax_inv,
                                     double stress_threshold_max_inv);

        hemelb::lb::LbmParameters *GetLbmParams();

        void RequestComms();
        void PreSend();
        void PreReceive();
        void PostReceive();
        void EndIteration();

        unsigned int siteMins[3], siteMaxes[3];

      private:
        void CalculateBC(double f[],
                         hemelb::geometry::LatticeData::SiteType iSiteType,
                         unsigned int iBoundaryId,
                         double *density,
                         double *vx,
                         double *vy,
                         double *vz,
                         double f_neq[]);

        void InitCollisions();

        //  static void ReadBlock();

        void ReadParameters();

        void allocateInlets(int nInlets);
        void allocateOutlets(int nOutlets);

        void handleIOError(int iError);

        double ConvertPressureGradToLatticeUnits(double pressure_grad) const;
        double ConvertPressureGradToPhysicalUnits(double pressure_grad) const;

        hemelb::lb::collisions::MidFluidCollision* mMidFluidCollision;
        hemelb::lb::collisions::WallCollision* mWallCollision;
        hemelb::lb::collisions::InletOutletCollision* mInletCollision;
        hemelb::lb::collisions::InletOutletCollision* mOutletCollision;
        hemelb::lb::collisions::InletOutletWallCollision* mInletWallCollision;
        hemelb::lb::collisions::InletOutletWallCollision* mOutletWallCollision;

        //TODO Get rid of this hack
        hemelb::lb::collisions::Collision* GetCollision(int i);

        double *inlet_density_avg, *inlet_density_amp;
        double *outlet_density_avg, *outlet_density_amp;
        double *inlet_density_phs, *outlet_density_phs;
        int is_inlet_normal_available;
        double* inlet_density, *outlet_density;
        hemelb::lb::collisions::MinsAndMaxes mMinsAndMaxes;
        double *inlet_normal;
        long int *inlet_count;
        double voxel_size;

        double mFileReadTime;

        SimConfig *mSimConfig;
        net::Net* mNet;
        geometry::LatticeData* mLatDat;
        SimulationState* mState;
        const topology::NetworkTopology * mNetTopology;

        LbmParameters mParams;
        vis::Control* mVisControl;

        double *average_inlet_velocity;
        double *peak_inlet_velocity;

        int * receivedFTranslator;
    };
  }
}
#endif // HEMELB_LB_H
