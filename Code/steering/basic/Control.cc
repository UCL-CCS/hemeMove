/*
 * Control.cc
 *
 *  Created on: Oct 27, 2010
 *      Author: rupert
 */
#include <pthread.h>
#include <semaphore.h>

#include "lb/lb.h"
#include "steering/Control.h"
#include "steering/basic/ClientConnection.h"
#include "steering/common/common.h"

namespace hemelb
{

  namespace steering
  {

    Control::Control(bool isCurrentProcTheSteeringProc) :
      is_frame_ready(false), sending_frame(false), isConnected(false), updated_mouse_coords(false),
          mIsCurrentProcTheSteeringProc(isCurrentProcTheSteeringProc)
    {
      // From main()
      sem_init(&nrl, 0, 1);
      sem_init(&steering_var_lock, 0, 1);
    }

    Control::~Control()
    {
      if (mIsCurrentProcTheSteeringProc)
      {
        delete mNetworkThread;
      }
      sem_destroy(&nrl);
      sem_destroy(&steering_var_lock);
    }

    // Kick off the networking thread
    void Control::StartNetworkThread(lb::LBM* lbm,
                                     lb::SimulationState *iSimState,
                                     const lb::LbmParameters *iLbmParams)
    {
      if (mIsCurrentProcTheSteeringProc)
      {
        ClientConnection* clientConnection = new ClientConnection(lbm->steering_session_id);
        mNetworkThread = new NetworkThread(lbm, this, iSimState, iLbmParams, clientConnection);
        mNetworkThread->Run();
      }
    }

    // Broadcast the steerable parameters to all tasks.
    void Control::UpdateSteerableParameters(bool shouldRenderForSnapshot,
                                            hemelb::lb::SimulationState &iSimulationState,
                                            hemelb::vis::Control* visController,
                                            lb::LBM* lbm)
    {
      if (mIsCurrentProcTheSteeringProc)
      {
        iSimulationState.DoRendering = ShouldRenderForNetwork() || shouldRenderForSnapshot;
        sem_wait(&steering_var_lock);
      }

      BroadcastSteerableParameters(iSimulationState, visController, lbm);

      if (mIsCurrentProcTheSteeringProc)
        sem_post(&steering_var_lock);
    }

    // Broadcast the steerable parameters to all tasks.
    void Control::BroadcastSteerableParameters(hemelb::lb::SimulationState &lSimulationState,
                                               vis::Control *visControl,
                                               lb::LBM* lbm)
    {
      steer_par[STEERABLE_PARAMETERS] = lSimulationState.DoRendering;

      MPI_Bcast(steer_par, STEERABLE_PARAMETERS + 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

      float longitude, latitude;
      float zoom;
      float lattice_velocity_max, lattice_stress_max;
      float lattice_density_min, lattice_density_max;

      int pixels_x, pixels_y;

      visControl->ctr_x += steer_par[0];
      visControl->ctr_y += steer_par[1];
      visControl->ctr_z += steer_par[2];

      longitude = steer_par[3];
      latitude = steer_par[4];

      zoom = steer_par[5];

      visControl->brightness = steer_par[6];

      // The minimum value here is by default 0.0 all the time
      visControl->physical_velocity_threshold_max = steer_par[7];

      // The minimum value here is by default 0.0 all the time
      visControl->physical_stress_threshold_max = steer_par[8];

      visControl->physical_pressure_threshold_min = steer_par[9];
      visControl->physical_pressure_threshold_max = steer_par[10];

      vis::GlyphDrawer::glyph_length = steer_par[11];

      pixels_x = steer_par[12];
      pixels_y = steer_par[13];

      int newMouseX = int (steer_par[14]);
      int newMouseY = int (steer_par[15]);

      if (newMouseX != visControl->mouse_x || newMouseY != visControl->mouse_y)
      {
        updated_mouse_coords = true;
        visControl->mouse_x = newMouseX;
        visControl->mouse_y = newMouseY;
      }

      lSimulationState.IsTerminating = int (steer_par[16]);

      // To swap between glyphs and streak line rendering...
      // 0 - Only display the isosurfaces (wall pressure and stress)
      // 1 - Isosurface and glyphs
      // 2 - Wall pattern streak lines
      visControl->mode = int (steer_par[17]);

      visControl->streaklines_per_pulsatile_period = steer_par[18];
      visControl->streakline_length = steer_par[19];

      lSimulationState.DoRendering = int (steer_par[20]);

      visControl->updateImageSize(pixels_x, pixels_y);

      lattice_density_min
          = lbm->ConvertPressureToLatticeUnits(visControl->physical_pressure_threshold_min) / Cs2;
      lattice_density_max
          = lbm->ConvertPressureToLatticeUnits(visControl->physical_pressure_threshold_max) / Cs2;
      lattice_velocity_max
          = lbm->ConvertVelocityToLatticeUnits(visControl->physical_velocity_threshold_max);
      lattice_stress_max
          = lbm->ConvertStressToLatticeUnits(visControl->physical_stress_threshold_max);

      visControl->SetProjection(pixels_x, pixels_y, visControl->ctr_x, visControl->ctr_y,
                                visControl->ctr_z, longitude, latitude, zoom);

      visControl->density_threshold_min = lattice_density_min;
      visControl->density_threshold_minmax_inv = 1.0F / (lattice_density_max - lattice_density_min);
      visControl->velocity_threshold_max_inv = 1.0F / lattice_velocity_max;
      visControl->stress_threshold_max_inv = 1.0F / lattice_stress_max;
    }

    // Do we need to render a frame for the client?
    bool Control::ShouldRenderForNetwork()
    {
      // Iff we're connected and not sending a
      // frame, we'd better render a new one!
      return isConnected.GetValue() && !sending_frame;
    }

    bool Control::RequiresSeparateSteeringCore() const
    {
      return true;
    }
  }

}
