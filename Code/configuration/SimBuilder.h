// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_CONFIGURATION_SIMBUILDER_H
#define HEMELB_CONFIGURATION_SIMBUILDER_H

#include "SimulationController.h"
#include "configuration/SimConfig.h"
#include "extraction/LbDataSourceIterator.h"
#include "extraction/PropertyActor.h"
#include "io/Checkpointer.h"
#include "geometry/GmyReadResult.h"
#include "geometry/neighbouring/NeighbouringDataManager.h"
#include "lb/lb.hpp"
#include "lb/InitialCondition.h"
#include "lb/StabilityTester.h"
#include "lb/IncompressibilityChecker.hpp"
#include "lb/iolets/BoundaryValues.h"
#include "net/PhasedBroadcastRegular.h"
#include "net/phased/StepManager.h"
#include "net/phased/NetConcern.h"
#include "reporting/Timers.h"
#include "util/clone_ptr.h"
#include "util/UnitConverter.h"

#ifdef HEMELB_BUILD_RBC
#include "redblood/CellControllerBuilder.h"
#endif

namespace hemelb::extraction { class PropertyActor; }
namespace hemelb::io { class Checkpointer; }
namespace hemelb::lb {
    class InOutLet;
}
namespace hemelb::net {
    class IteratedAction;
}
namespace hemelb::reporting {
    class BuildInfo;
    class Reporter;
}

namespace hemelb::configuration {

    // This class converters the SimConfig into a ready-to-run simulation.
    class SimBuilder {
    public:
        using IoletPtr = util::clone_ptr<lb::InOutLet>;

    protected:
        SimConfig config;
        std::shared_ptr<util::UnitConverter> unit_converter;

    public:
        template <typename TraitsT>
        static std::unique_ptr<SimulationController> CreateSim(
                CommandLine const& options,
                net::IOCommunicator const& ioComms
        ) {
            auto ans = std::unique_ptr<SimulationController>(new SimulationController(ioComms));
            // Start the main timer!
            ans->timings.total().Start();

            ans->fileManager = std::make_shared<configuration::PathManager>(
                    options,
                    ioComms.OnIORank(),
                    ioComms.Size()
            );
            auto&& infile = ans->fileManager->GetInputFile();
            log::Logger::Log<log::Info, log::Singleton>("Reading configuration from %s", infile.c_str());
            // Convert XML to configuration
            ans->simConfig = configuration::SimConfig::New(infile);
            // Use it to initialise self
            auto builder = configuration::SimBuilder(ans->simConfig);
            log::Logger::Log<log::Info, log::Singleton>("Beginning Initialisation.");
            builder.build<TraitsT>(*ans);
            return ans;
        }

        explicit SimBuilder(SimConfig const& conf, bool construct_unit_converter = true);
        virtual ~SimBuilder() = default;

        [[nodiscard]] std::shared_ptr<util::UnitConverter const> GetUnitConverter() const;

        template<typename T>
        T ConvertToLatticeUnits(T const& val, std::string_view units)
        {
            return unit_converter->template ConvertToLatticeUnits(units, val);
        }

        // Fully build the T = SimulationController<Traits> from the configuration.
        template <typename TRAITS = hemelb::Traits<>>
        void build(SimulationController& control) const;

        // The below could probably be protected/private, but handy for testing.
        [[nodiscard]] std::shared_ptr<lb::SimulationState> BuildSimulationState() const;
        [[nodiscard]] geometry::GmyReadResult ReadGmy(
                lb::LatticeInfo const& lat_info,
                reporting::Timers& timings,
                net::IOCommunicator& ioComms
        ) const;

        [[nodiscard]] lb::LbmParameters BuildLbmParams() const;

        [[nodiscard]] std::vector<IoletPtr> BuildIolets(std::vector<IoletConfig> const&) const;
        [[nodiscard]] IoletPtr BuildIolet(IoletConfig const&) const;
        [[nodiscard]] IoletPtr BuildCosinePressureIolet(CosinePressureIoletConfig const&) const;
        [[nodiscard]] IoletPtr BuildFilePressureIolet(FilePressureIoletConfig const&) const;
        [[nodiscard]] IoletPtr BuildMultiscalePressureIolet(MultiscalePressureIoletConfig const&) const;
        [[nodiscard]] IoletPtr BuildParabolicVelocityIolet(ParabolicVelocityIoletConfig const&) const;
        [[nodiscard]] IoletPtr BuildWomersleyVelocityIolet(WomersleyVelocityIoletConfig const&) const;
        [[nodiscard]] IoletPtr BuildFileVelocityIolet(FileVelocityIoletConfig const&) const;
        void BuildBaseIolet(IoletConfigBase const& conf, lb::InOutLet* obj) const;
        [[nodiscard]] std::shared_ptr<redblood::FlowExtension> BuildFlowExtension(FlowExtensionConfig const& conf) const;

        [[nodiscard]] std::shared_ptr<net::IteratedAction> BuildColloidController() const;
        template <typename traitsType>
        [[nodiscard]] std::shared_ptr<net::IteratedAction> BuildCellController(SimulationController const& control, reporting::Timers& timings) const;
        [[nodiscard]] lb::InitialCondition BuildInitialCondition() const;
        [[nodiscard]] std::shared_ptr<extraction::PropertyActor> BuildPropertyExtraction(
                std::filesystem::path const& xtr_path,
                std::shared_ptr<lb::SimulationState const> simState,
                std::shared_ptr<extraction::IterableDataSource> dataSource,
                reporting::Timers& timings,
                const net::IOCommunicator& ioComms
        ) const;

        [[nodiscard]] std::shared_ptr<io::Checkpointer> BuildCheckpointer(
                std::filesystem::path const& cp_path,
                std::shared_ptr<lb::SimulationState const> simState,
                std::shared_ptr<extraction::IterableDataSource> dataSource,
                reporting::Timers& timings,
                const net::IOCommunicator& ioComms
        ) const;

        [[nodiscard]] std::shared_ptr<reporting::Reporter> BuildReporter(
                PathManager const& fileManager,
                std::vector<reporting::Reportable*> const & reps
        ) const;
    };


    template <typename traitsType>
    void SimBuilder::build(SimulationController& control) const {
        using LatticeType = typename traitsType::Lattice;

        auto& timings = control.timings;
        auto& ioComms = control.ioComms;
        auto& lat_info = LatticeType::GetLatticeInfo();

        control.unitConverter = unit_converter;

        control.simulationState = BuildSimulationState();

        std::vector<reporting::Reportable*> things_to_report({
            &control.build_info, &timings, &*control.simulationState
        });

        std::vector<std::pair<net::phased::Concern*, unsigned>> actors_to_register_for_phase;
        auto maybe_register_actor = [&] <std::derived_from<net::phased::Concern> C> (std::shared_ptr<C> const& p, unsigned i) {
            if (p)
                actors_to_register_for_phase.emplace_back(p.get(), i);
        };

        timings.latDatInitialise().Start();
        // Use a reader to read in the file.
        log::Logger::Log<log::Info, log::Singleton>("Loading and decomposing geometry file %s.", config.GetDataFilePath().c_str());
        auto readGeometryData = ReadGmy(lat_info, timings, ioComms);
        // Create a new lattice based on that info and return it.
        log::Logger::Log<log::Info, log::Singleton>("Initialising domain.");
        control.domainData = std::make_shared<geometry::Domain>(lat_info,
                                                                readGeometryData,
                                                                ioComms);
        log::Logger::Log<log::Info, log::Singleton>("Initialising field data.");
        control.fieldData = std::make_shared<geometry::FieldData>(control.domainData);
        things_to_report.push_back(control.domainData.get());
        timings.latDatInitialise().Stop();

        log::Logger::Log<log::Info, log::Singleton>("Initialising neighbouring data manager.");
        auto ndm =
                control.neighbouringDataManager =
                        std::make_shared<geometry::neighbouring::NeighbouringDataManager>(
                                *control.fieldData,
                                control.fieldData->GetNeighbouringData(),
                                control.communicationNet
                        );
        maybe_register_actor(ndm, 0);

        log::Logger::Log<log::Info, log::Singleton>("Initialising LBM.");
        auto lbm =
                control.latticeBoltzmannModel =
                        std::make_shared<lb::LBM<traitsType>>(
                                BuildLbmParams(),
                                &control.communicationNet,
                                control.fieldData.get(),
                                &*control.simulationState,
                                timings,
                                control.neighbouringDataManager.get()
                        );

        maybe_register_actor(control.colloidController = BuildColloidController(), 1);

        maybe_register_actor(lbm, 1);

        control.inletValues = std::make_shared<lb::BoundaryValues>(
                geometry::INLET_TYPE,
                *control.domainData,
                BuildIolets(config.GetInlets()),
                &*control.simulationState,
                ioComms,
                *unit_converter
        );
        maybe_register_actor(control.inletValues, 1);

        control.outletValues = std::make_shared<lb::BoundaryValues>(
                geometry::OUTLET_TYPE,
                *control.domainData,
                BuildIolets(config.GetOutlets()),
                &*control.simulationState,
                ioComms,
                *unit_converter
        );
        maybe_register_actor(control.outletValues, 1);

        maybe_register_actor(control.cellController = BuildCellController<traitsType>(control, timings), 1);

        // Copy cos about to scale to lattice units.
        auto mon_conf = config.GetMonitoringConfiguration();
        if (std::holds_alternative<extraction::source::Velocity>(mon_conf.convergenceVariable)) {
            mon_conf.convergenceReferenceValue = unit_converter->ConvertSpeedToLatticeUnits(mon_conf.convergenceReferenceValue);
        }

        // Always track stability
        control.stabilityTester = std::make_shared<lb::StabilityTesterImpl<LatticeType>>(
                control.fieldData,
                &control.communicationNet,
                &*control.simulationState,
                timings,
                mon_conf
        );
        maybe_register_actor(control.stabilityTester, 1);

        // Incompressibility only if requested
        if (mon_conf.doIncompressibilityCheck)
        {
            control.incompressibilityChecker =
                    std::make_shared<lb::IncompressibilityChecker<net::PhasedBroadcastRegular<>>>(
                            control.domainData.get(),
                            &control.communicationNet,
                            &*control.simulationState,
                            control.latticeBoltzmannModel->GetPropertyCache(),
                            timings
                    );
            things_to_report.push_back(control.incompressibilityChecker.get());
            maybe_register_actor(control.incompressibilityChecker, 1);
        }

        lbm->Initialise(control.inletValues.get(),
                        control.outletValues.get());
        auto ic = BuildInitialCondition();
        lbm->SetInitialConditions(ic, ioComms);
        ndm->ShareNeeds();
        ndm->TransferNonFieldDependentInformation();

        control.propertyDataSource = std::make_shared<extraction::LbDataSourceIterator>(
                lbm->GetPropertyCache(),
                *control.fieldData,
                ioComms.Rank(),
                unit_converter
        );

        control.propertyExtractor = BuildPropertyExtraction(
                control.fileManager->GetDataExtractionPath(),
                control.simulationState,
                control.propertyDataSource,
                timings,
                ioComms
        );
        maybe_register_actor(control.propertyExtractor, 1);

        control.checkpointer = BuildCheckpointer(
                control.fileManager->GetCheckpointPath(),
                control.simulationState,
                control.propertyDataSource,
                timings,
                ioComms
        );

        control.netConcern = std::make_shared<net::phased::NetConcern>(
                control.communicationNet
        );

        control.reporter = BuildReporter(*control.fileManager, things_to_report);

        control.stepManager = std::make_shared<net::phased::StepManager>(
                2,
                &timings,
                net::separate_communications
        );
        for (auto [actor_ptr, phase]: actors_to_register_for_phase) {
            control.stepManager->RegisterIteratedActorSteps(*actor_ptr, phase);
        }
        control.stepManager->RegisterCommsForAllPhases(*control.netConcern);
    }

#ifdef HEMELB_BUILD_RBC
    inline redblood::CountedIoletView MakeCountedIoletView(lb::BoundaryValues const& iolets) {
        return {
                [&]() { return iolets.GetGlobalIoletCount(); },
                [&](unsigned i) { return iolets.GetGlobalIolet(i); }
        };
    }
    inline redblood::CountedIoletView MakeCountedIoletView(std::vector<SimBuilder::IoletPtr> const& iolets) {
        return {
                [&]() { return iolets.size(); },
                [&](unsigned i) { return iolets[i].get(); }
        };
    }
#endif

    template <typename traitsType>
    [[nodiscard]] std::shared_ptr<net::IteratedAction> SimBuilder::BuildCellController(SimulationController const& control, reporting::Timers& timings) const {
        if (config.HasRBCSection()) {
#ifdef HEMELB_BUILD_RBC
            auto ccb = redblood::CellControllerBuilder(unit_converter);

            auto& ioComms = control.ioComms;

            auto inlets = MakeCountedIoletView(*control.inletValues);
            auto outlets = MakeCountedIoletView(*control.outletValues);

            return ccb.build<traitsType>(
                    config,
                    timings, ioComms,
                    control.fieldData,
                    inlets,
                    outlets,
                    control.simulationState,
                    control.fileManager
            );
#else
            throw hemelb::Exception() << "Trying to create red blood cell controller with HEMELB_BUILD_RBC=OFF";
#endif
        }
        return {};
    }
}
#endif
