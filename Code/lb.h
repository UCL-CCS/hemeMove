#ifndef HEMELB_LB_H
#define HEMELB_LB_H

#include "net.h"
#include "lb/collisions/Collisions.h"
#include "vis/ColPixel.h"
#include "SimConfig.h"

class LBM
{
  public:
    double *inlet_density_avg, *inlet_density_amp;
    double *outlet_density_avg, *outlet_density_amp;

    int total_fluid_sites;
    int inlets, outlets;
    int steering_session_id;
    double voxel_size;
    int period;

    double lbmConvertPressureToLatticeUnits(double pressure);
    double lbmConvertPressureToPhysicalUnits(double pressure);
    double lbmConvertVelocityToLatticeUnits(double velocity);
    double lbmConvertStressToLatticeUnits(double stress);
    double lbmConvertStressToPhysicalUnits(double stress);
    double lbmConvertVelocityToPhysicalUnits(double velocity);

    void lbmInit(hemelb::SimConfig *iSimulationConfig,
                 int iSteeringSessionId,
                 int iPeriod,
                 double voxel_size,
                 Net *net);
    void lbmRestart(Net *net);
    ~LBM();

    int IsUnstable(Net *net);

    int lbmCycle(int perform_rt, Net *net);
    void lbmCalculateFlowFieldValues();
    void RecalculateTauViscosityOmega();
    void lbmUpdateBoundaryDensities(int cycle_id, int time_step);
    void lbmUpdateInletVelocities(int time_step, Net *net);

    void lbmSetInitialConditions(Net *net);

    void lbmWriteConfig(int stability, std::string output_file_name, Net *net);

    double GetMinPhysicalPressure();
    double GetMaxPhysicalPressure();
    double GetMinPhysicalVelocity();
    double GetMaxPhysicalVelocity();
    double GetMinPhysicalStress();
    double GetMaxPhysicalStress();

    void lbmInitMinMaxValues(void);

    double GetAverageInletVelocity(int iInletNumber);
    double GetPeakInletVelocity(int iInletNumber);

    void ReadVisParameters(Net *net, hemelb::vis::Control *bController);

    void CalculateMouseFlowField(hemelb::vis::ColPixel *col_pixel_p,
                                 double &mouse_pressure,
                                 double &mouse_stress,
                                 double density_threshold_min,
                                 double density_threshold_minmax_inv,
                                 double stress_threshold_max_inv);

  private:
    void lbmCalculateBC(double f[],
                        unsigned int site_data,
                        double *density,
                        double *vx,
                        double *vy,
                        double *vz,
                        double f_neq[]);

    void lbmInitCollisions();
    void lbmReadConfig(Net *net);
    void lbmReadParameters(Net *net);

    void allocateInlets(int nInlets);
    void allocateOutlets(int nOutlets);

    double lbmConvertPressureGradToLatticeUnits(double pressure_grad);
    double lbmConvertPressureGradToPhysicalUnits(double pressure_grad);

    hemelb::lb::collisions::MidFluidCollision* mMidFluidCollision;
    hemelb::lb::collisions::WallCollision* mWallCollision;
    hemelb::lb::collisions::InletOutletCollision* mInletCollision;
    hemelb::lb::collisions::InletOutletCollision* mOutletCollision;
    hemelb::lb::collisions::InletOutletWallCollision* mInletWallCollision;
    hemelb::lb::collisions::InletOutletWallCollision* mOutletWallCollision;

    //TODO Get rid of this hack
    hemelb::lb::collisions::Collision* GetCollision(int i);

    double *inlet_density_phs, *outlet_density_phs;

    int site_min_x, site_min_y, site_min_z;
    int site_max_x, site_max_y, site_max_z;
    int is_inlet_normal_available;
    double* inlet_density, *outlet_density;
    hemelb::lb::collisions::MinsAndMaxes mMinsAndMaxes;
    double *lbm_inlet_normal;
    long int *lbm_inlet_count;
    double tau, viscosity, omega;

    // TODO Eventually we should be able to make this private.
    hemelb::SimConfig *mSimConfig;

    double *lbm_average_inlet_velocity;
    double *lbm_peak_inlet_velocity;
};

extern double lbm_stress_type;
extern double lbm_stress_par;
extern int lbm_terminate_simulation;

#endif // HEMELB_LB_H
