#ifndef HEMELB_LB_H
#define HEMELB_LB_H

#include "net/net.h"
#include "topology/NetworkTopology.h"
#include "lb/collisions/Collisions.h"
#include "vis/ColPixel.h"
#include "SimConfig.h"

namespace hemelb
{
  namespace lb
  {
    class LBM
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
                const hemelb::topology::NetworkTopology * iNetTop);
        void Restart(geometry::LatticeData* iLatDat);
        ~LBM();

        void CalculateFlowFieldValues();
        void RecalculateTauViscosityOmega();
        void UpdateBoundaryDensities(int cycle_id, int time_step);
        void
        UpdateInletVelocities(int time_step, geometry::LatticeData &iLatDat, net::Net *net);

        void Initialise(int* iFTranslator, geometry::LatticeData* bLatDat, vis::Control* iControl);

        void SetInitialConditions(geometry::LatticeData* bLatDat);

        void
        WriteConfig(hemelb::lb::Stability stability,
                    std::string output_file_name,
                    const geometry::LatticeData &iLatticeData);
        void
        WriteConfigParallel(hemelb::lb::Stability stability,
                            std::string output_file_name,
                            const geometry::LatticeData &iLatticeData);

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

        void RequestComms(net::Net* net, geometry::LatticeData* bLatDat);
        void PreSend(geometry::LatticeData* bLatDat, int perform_rt);
        void PreReceive(int perform_rt, geometry::LatticeData* bLatDat);
        void PostReceive(geometry::LatticeData* bLatDat, int perform_rt);
        void EndIteration(geometry::LatticeData* bLatDat);

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

        LbmParameters mParams;
        const topology::NetworkTopology * mNetTopology;
        SimConfig *mSimConfig;
        vis::Control* mVisControl;

        double *average_inlet_velocity;
        double *peak_inlet_velocity;

        int * receivedFTranslator;
    };
  }
}
#endif // HEMELB_LB_H
