// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_SIMULATIONCONTROLLER_H
#define HEMELB_SIMULATIONCONTROLLER_H
#include <memory>

#include "lb/Lattices.h"
#include "lb/lb.h"
#include "lb/StabilityTester.h"
#include "net/net.h"
#include "lb/EntropyTester.h"
#include "lb/iolets/BoundaryValues.h"
#include "util/UnitConverter.h"
#include "configuration/SimConfig.h"
#include "configuration/CommandLine.h"
#include "configuration/PathManager.h"
#include "reporting/Reporter.h"
#include "reporting/Timers.h"
#include "reporting/BuildInfo.h"
#include "lb/IncompressibilityChecker.hpp"
#include "net/phased/StepManager.h"
#include "net/phased/NetConcern.h"
#include "geometry/neighbouring/NeighbouringDataManager.h"
#include "Traits.h"

namespace hemelb
{
    namespace configuration { class SimBuilder; }
    namespace extraction { class PropertyActor; }
    namespace io { class Checkpointer; }

    class SimulationController
    {
    public:
        friend class configuration::SimBuilder;

    protected:
        reporting::BuildInfo build_info;
        reporting::Timers timings;
        net::IOCommunicator ioComms;
        net::Net communicationNet;

        configuration::SimConfig simConfig;

        std::shared_ptr<lb::BoundaryValues> inletValues;
        std::shared_ptr<lb::BoundaryValues> outletValues;
        /* The next quantities are protected because they are used by MultiscaleSimulationController */
        std::shared_ptr<geometry::Domain> domainData;
        std::shared_ptr<geometry::FieldData> fieldData;
        std::shared_ptr<lb::LBMBase> latticeBoltzmannModel;
        std::shared_ptr<geometry::neighbouring::NeighbouringDataManager> neighbouringDataManager;

        std::shared_ptr<configuration::PathManager> fileManager;
        std::shared_ptr<reporting::Reporter> reporter;

        std::shared_ptr<lb::SimulationState> simulationState;

        /** Struct containing the configuration of various checkers/testers */
        std::shared_ptr<lb::StabilityTester> stabilityTester;
        // std::shared_ptr<lb::EntropyTester<LatticeType>> entropyTester;
        /** Actor in charge of checking the maximum density difference across the domain */
        using ICC = lb::IncompressibilityChecker<net::PhasedBroadcastRegular<>>;
        std::shared_ptr<ICC> incompressibilityChecker;

        std::shared_ptr<net::IteratedAction> cellController;
        std::shared_ptr<net::IteratedAction> colloidController;

        std::shared_ptr<util::UnitConverter> unitConverter;

        std::shared_ptr<extraction::IterableDataSource> propertyDataSource;
        std::shared_ptr<extraction::PropertyActor> propertyExtractor;
        std::shared_ptr<io::Checkpointer> checkpointer;

        std::shared_ptr<net::phased::StepManager> stepManager;
        std::shared_ptr<net::phased::NetConcern> netConcern;

        SimulationController(const net::IOCommunicator& ioComm);
    public:
      virtual ~SimulationController();

      void Abort();

      /// Returns true if the current processor is the dedicated I/O processor.
      inline bool IsCurrentProcTheIOProc() {
          return ioComms.OnIORank();
      }

      int GetProcessorCount();

      void RunSimulation();
      lb::SimulationState const& GetState() const
      {
        return *simulationState;
      }
      util::UnitConverter const & GetUnitConverter() const;
      void Finalise();

      configuration::SimConfig ToConfig() const;
#     ifdef HEMELB_DOING_UNITTESTS
      //! Makes it easy to add cell controller without messy input files
      void RegisterActor(net::phased::Concern &concern, net::phased::StepManager::Phase phase)
      {
        stepManager->RegisterIteratedActorSteps(concern, phase);
      }
      //! Access to lattice data for debugging
      geometry::FieldData& GetFieldData()
      {
        assert(fieldData);
        return *fieldData;
      }
      std::shared_ptr<net::IteratedAction> GetCellController() {
        return cellController;
      }
#     endif
    protected:

      virtual void DoTimeStep();


    private:
      void Initialise();
      void SetupReporting(); // set up the reporting file
      void HandleActors();
      void OnUnstableSimulation();
      /**
       * Updates the property caches record of which properties need to be calculated
       * and cached on this iteration.
       */
      void RecalculatePropertyRequirements();

      /**
       * Helper method to log simulation parameters related to stability and accuracy
       */
      void LogStabilityReport();



      static constexpr LatticeTimeStep FORCE_FLUSH_PERIOD = 1000;
  };
}

#endif /* HEMELB_SIMULATIONCONTROLLER_H */
