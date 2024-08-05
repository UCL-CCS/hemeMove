// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "SimulationController.h"

#include <map>
#include <limits>
#include <cstdlib>
#include <boost/uuid/uuid_io.hpp>

#include "configuration/SimConfig.h"
#include "configuration/SimBuilder.h"
#include "extraction/PropertyActor.h"
#include "extraction/LbDataSourceIterator.h"
#include "io/writers/XdrFileWriter.h"
#include "util/numerical.h"
#include "geometry/Domain.h"
#include "log/Logger.h"
#include "lb/HFunction.h"
#include "io/xml.h"
#include "net/BuildInfo.h"
#include "net/IOCommunicator.h"

#ifdef HEMELB_BUILD_COLLOIDS
#  include "colloids/BodyForces.h"
#  include "colloids/BoundaryConditions.h"
#  include "colloids/ColloidController.h"
#endif

namespace hemelb
{
    /**
     * Constructor for the SimulationController class
     *
     * Initialises member variables including the network topology
     * object.
     */
    SimulationController::SimulationController(const net::IOCommunicator& ioComm) :
            build_info(), ioComms(ioComm.Duplicate()), communicationNet(ioComms)
    {
    }

  /// Destructor for the SimulationController class.
  SimulationController::~SimulationController() = default;

  /**
   * Returns the number of processors involved in the simulation.
   */
  int SimulationController::GetProcessorCount()
  {
    return ioComms.Size();
  }

  /**
   * Initialises various elements of the simulation if necessary:
   * domain decomposition, LBM and visualisation.
   */
  void SimulationController::Initialise()
  {

  }

  void SimulationController::HandleActors()
  {
    stepManager->CallActions();
  }

  void SimulationController::OnUnstableSimulation()
  {
    LogStabilityReport();
    log::Logger::Log<log::Warning, log::Singleton>("Aborting: time step length: %f\n",
                                                                           simulationState->GetTimeStepLength());
    Finalise();
    Abort();
  }

  /**
   * Begin the simulation.
   */
  void SimulationController::RunSimulation()
  {
    log::Logger::Log<log::Info, log::Singleton>("Beginning to run simulation.");
    timings.simulation().Start();

    auto cp = [&] () {
        if (checkpointer && checkpointer->ShouldWrite()) {
            timings.writeCheckpoint().Start();
            checkpointer->Write(ToConfig());
            timings.writeCheckpoint().Stop();
        }
    };

    cp();
    while (simulationState->GetTimeStep() < simulationState->GetEndTimeStep())
    {
      DoTimeStep();
      cp();

      if (simulationState->IsTerminating())
      {
        break;
      }
    }

    timings.simulation().Stop();
    Finalise();
  }

  void SimulationController::Finalise()
  {
    timings.total().Stop();
    timings.Reduce(ioComms);

    if (IsCurrentProcTheIOProc())
    {
      reporter->FillDictionary();
      reporter->Write();
    }
    // DTMP: Logging output on communication as debug output for now.
    log::Logger::Log<log::Debug, log::OnePerCore>("sync points: %lld, bytes sent: %lld",
                                                                          communicationNet.SyncPointsCounted,
                                                                          communicationNet.BytesSent);

    log::Logger::Log<log::Info, log::Singleton>("Finish running simulation.");
  }

  void SimulationController::DoTimeStep()
  {
    log::Logger::Log<log::Debug, log::OnePerCore>("Current LB time: %e",
                                                  simulationState->GetTime());

    /* In the following two if blocks we do the core magic to ensure we only Render
     when (1) we are not sending a frame or (2) we need to output to disk */

    /* TODO for debugging purposes we want to ensure we capture the variables in a single
     instant of time since variables might be altered by the thread half way through?
     This is to be done. */

    if (simulationState->GetTimeStep() % 100 == 0)
    {
      log::Logger::Log<log::Info, log::Singleton>("time step %i",
                                                                          simulationState->GetTimeStep());
      LogStabilityReport();
    }

    RecalculatePropertyRequirements();

    HandleActors();

    if (simulationState->GetStability() == lb::Unstable)
    {
      OnUnstableSimulation();
    }

    // If the user requested to terminate converged steady flow simulations, mark
    // simulation to be finished at the end of the current timestep.
    if ( (simulationState->GetStability() == lb::StableAndConverged)
        && stabilityTester->ShouldTerminateWhenConverged())
    {
      LogStabilityReport();
      simulationState->SetIsTerminating(true);
    }

#ifdef HEMELB_BUILD_COLLOIDS
    if (simulationState->GetTimeStep() % 500 == 0) {
        if (auto c = std::dynamic_pointer_cast<colloids::ColloidController>(colloidController)) {
            c->OutputInformation(simulationState->GetTimeStep());
        }
    }
#endif
    if (simulationState->GetTimeStep() % FORCE_FLUSH_PERIOD == 0 && IsCurrentProcTheIOProc())
    {
      fflush(nullptr);
    }

    fieldData->SwapOldAndNew();
    simulationState->Increment();
  }

  void SimulationController::RecalculatePropertyRequirements()
  {
    // Get the property cache & reset its list of properties to get.
    lb::MacroscopicPropertyCache& propertyCache = latticeBoltzmannModel->GetPropertyCache();

    propertyCache.ResetRequirements();

    if (incompressibilityChecker)
    {
      propertyCache.densityCache.SetRefreshFlag();
      propertyCache.velocityCache.SetRefreshFlag();
    }

    // If extracting property results, check what's required by them.
    if (propertyExtractor)
    {
      propertyExtractor->SetRequiredProperties(propertyCache);
    }
  }

  /**
   * Called on error to abort the simulation and pull-down the MPI environment.
   */
  void SimulationController::Abort()
  {
    // This gives us something to work from when we have an error - we get the rank
    // that calls abort, and we get a stack-trace from the exception having been thrown.
    log::Logger::Log<log::Critical, log::OnePerCore>("Aborting");
    net::MpiEnvironment::Abort(1);

    exit(1);
  }

  void SimulationController::LogStabilityReport()
  {
    if (incompressibilityChecker && incompressibilityChecker->AreDensitiesAvailable())
    {
      log::Logger::Log<log::Info, log::Singleton>("time step %i, tau %.6f, max_relative_press_diff %.3f, Ma %.3f, max_vel_phys %e",
                                                                          simulationState->GetTimeStep(),
                                                                          latticeBoltzmannModel->GetLbmParams()->GetTau(),
                                                                          incompressibilityChecker->GetMaxRelativeDensityDifference(),
                                                                          incompressibilityChecker->GetGlobalLargestVelocityMagnitude()
                                                                              / Cs,
                                                                          unitConverter->ConvertVelocityToPhysicalUnits(incompressibilityChecker->GetGlobalLargestVelocityMagnitude()));
    }

    if (simulationState->GetStability() == lb::StableAndConverged)
    {
      log::Logger::Log<log::Info, log::Singleton>("time step %i, steady flow simulation converged.",
                                                                          simulationState->GetTimeStep());
    }
  }

  const util::UnitConverter& SimulationController::GetUnitConverter() const
  {
    return *unitConverter;
  }

    configuration::SimConfig SimulationController::ToConfig() const {
        using namespace configuration;
        SimConfig ans = simConfig;
        if (simConfig.HasColloidSection())
            throw (Exception() << "Checkpointing not implemented for colloids");
        return ans;
    }
}
