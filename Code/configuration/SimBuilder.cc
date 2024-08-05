// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "configuration/SimBuilder.h"

#include <algorithm>
#include <ranges>

#include "extraction/GeometrySelectors.h"
#include "geometry/GeometryReader.h"
#include "lb/InitialCondition.h"
#include "lb/iolets/InOutLets.h"
#include "redblood/FlowExtension.h"
#include "reporting/Reporter.h"
#include "util/variant.h"

namespace hemelb::configuration {

    auto build_unit_converter(GlobalSimInfo const& info) {
        return std::make_shared<util::UnitConverter>(
            info.time.step_s,
            info.space.step_m, info.space.geometry_origin_m,
            info.fluid.density_kgm3, info.fluid.reference_pressure_Pa
        );
    }

    SimBuilder::SimBuilder(configuration::SimConfig const& conf, bool construct_unit_converter) :
            config(conf)
    {
        if (construct_unit_converter)
            unit_converter = build_unit_converter(conf.sim_info);
    }

    std::shared_ptr<util::UnitConverter const> SimBuilder::GetUnitConverter() const {
        return unit_converter;
    }

    std::shared_ptr<lb::SimulationState> SimBuilder::BuildSimulationState() const {
        return std::make_shared<lb::SimulationState>(config.GetTimeStepLength(), config.GetTotalTimeSteps());
    }

    geometry::GmyReadResult SimBuilder::ReadGmy(lb::LatticeInfo const& lat_info, reporting::Timers& timings, net::IOCommunicator& ioComms) const {
        geometry::GeometryReader reader(lat_info,
                                        timings,
                                        ioComms);
        return reader.LoadAndDecompose(config.GetDataFilePath());
    }

    lb::LbmParameters SimBuilder::BuildLbmParams() const {
        auto&& i = config.sim_info;
        lb::LbmParameters ans(i.time.step_s, i.space.step_m, i.fluid.density_kgm3, i.fluid.viscosity_Pas);
        return ans;
    }

    // Visitor for factory function
    struct ICMaker {
        using result_type = lb::InitialCondition;
        util::UnitConverter const& units;

        template <typename T>
        result_type operator()(T) const {
            throw Exception() << "Trying to make an InitialCondition from unknown type of config";
        }

        result_type operator()(const configuration::EquilibriumIC& cfg) const {
            auto rho = units.ConvertPressureToLatticeUnits(cfg.p_Pa) / Cs2;
            return lb::EquilibriumInitialCondition{cfg.t0, rho};
        }
        result_type operator()(const configuration::CheckpointIC& cfg) const {
            return lb::CheckpointInitialCondition{cfg.t0, cfg.cpFile, cfg.maybeOffFile};
        }
    };

    // Factory function just delegates to visitor
    lb::InitialCondition SimBuilder::BuildInitialCondition() const {
        return std::visit(ICMaker{*unit_converter}, config.initial_condition);
    }

    // Build the iolets
    auto SimBuilder::BuildIolets(const std::vector<IoletConfig>& io_confs) const -> std::vector<IoletPtr> {
        std::vector<IoletPtr> ans;
        std::transform(io_confs.begin(), io_confs.end(), std::back_inserter(ans),
                       [&](IoletConfig const& ic) {
            return BuildIolet(ic);
        });
        return ans;
    }

    auto SimBuilder::BuildIolet(const IoletConfig& ic) const -> IoletPtr {
        return overload_visit(
                ic,
                [] (std::monostate) -> IoletPtr { throw Exception() << "Invalid IoletConfig"; },
                [&](CosinePressureIoletConfig const& _) { return BuildCosinePressureIolet(_); },
                [&](FilePressureIoletConfig const& _) { return BuildFilePressureIolet(_); },
                [&](MultiscalePressureIoletConfig const& _) { return BuildMultiscalePressureIolet(_); },
                [&](ParabolicVelocityIoletConfig const& _) { return BuildParabolicVelocityIolet(_); },
                [&](WomersleyVelocityIoletConfig const& _) { return BuildWomersleyVelocityIolet(_); },
                [&](FileVelocityIoletConfig const& _) { return BuildFileVelocityIolet(_); }
        );
    }
    auto SimBuilder::BuildCosinePressureIolet(const CosinePressureIoletConfig& ic) const -> IoletPtr {
        auto ans = util::make_clone_ptr<lb::InOutLetCosine>();
        BuildBaseIolet(ic, ans.get());

        // Amplitude is a pressure DIFFERENCE (no use of REFERENCE_PRESSURE)
        ans->SetPressureAmp(unit_converter->ConvertPressureDifferenceToLatticeUnits(ic.amp_Pa));
        // Mean is an absolute pressure
        ans->SetPressureMean(unit_converter->ConvertPressureToLatticeUnits(ic.mean_Pa));
        ans->SetPhase(ic.phase_rad);
        ans->SetPeriod(unit_converter->ConvertTimeToLatticeUnits(ic.period_s));
        return ans;
    }

    auto SimBuilder::BuildFilePressureIolet(const FilePressureIoletConfig & ic) const -> IoletPtr {
        auto ans = util::make_clone_ptr<lb::InOutLetFile>();
        BuildBaseIolet(ic, ans.get());
        ans->SetFilePath(ic.file_path);
        return ans;
    }

    auto SimBuilder::BuildMultiscalePressureIolet(const MultiscalePressureIoletConfig & ic) const -> IoletPtr {
        auto ans = util::make_clone_ptr<lb::InOutLetMultiscale>();
        BuildBaseIolet(ic, ans.get());
        ans->GetPressureReference() = ic.pressure_reference_Pa;
        ans->GetVelocityReference() = ic.velocity_reference_ms;
        ans->GetLabel() = ic.label;
        return ans;
    }

    auto SimBuilder::BuildParabolicVelocityIolet(const ParabolicVelocityIoletConfig& ic) const -> IoletPtr {
        auto ans = util::make_clone_ptr<lb::InOutLetParabolicVelocity>();
        BuildBaseIolet(ic, ans.get());
        ans->SetRadius(unit_converter->ConvertDistanceToLatticeUnits(ic.radius_m));
        ans->SetMaxSpeed(unit_converter->ConvertSpeedToLatticeUnits(ic.max_speed_ms));
        return ans;
    }

    auto SimBuilder::BuildWomersleyVelocityIolet(const WomersleyVelocityIoletConfig& ic) const -> IoletPtr {
        auto ans = util::make_clone_ptr<lb::InOutLetWomersleyVelocity>();
        BuildBaseIolet(ic, ans.get());
        ans->SetRadius(unit_converter->ConvertDistanceToLatticeUnits(ic.radius_m));
        ans->SetPressureGradientAmplitude(unit_converter->ConvertPressureGradientToLatticeUnits(ic.pgrad_amp_Pam));
        ans->SetPeriod(unit_converter->ConvertTimeToLatticeUnits(ic.period_s));
        ans->SetWomersleyNumber(ic.womersley);
        return ans;
    }

    auto SimBuilder::BuildFileVelocityIolet(const FileVelocityIoletConfig &ic) const -> IoletPtr {
        auto ans = util::make_clone_ptr<lb::InOutLetFileVelocity>();
        BuildBaseIolet(ic, ans.get());
        ans->SetFilePath(ic.file_path);
        ans->SetRadius(unit_converter->ConvertDistanceToLatticeUnits(ic.radius_m));
        return ans;
    }

    void SimBuilder::BuildBaseIolet(IoletConfigBase const& conf,
                                        lb::InOutLet* obj) const
    {
        obj->SetPosition(unit_converter->ConvertPositionToLatticeUnits(conf.position));
        obj->SetNormal(conf.normal);
        if (conf.flow_extension.has_value()) {
            obj->SetFlowExtension(BuildFlowExtension(conf.flow_extension.value()));
        }
    }

    std::shared_ptr<redblood::FlowExtension> SimBuilder::BuildFlowExtension(FlowExtensionConfig const& conf) const {
        return std::make_shared<redblood::FlowExtension>(
                conf.normal,
                unit_converter->ConvertPositionToLatticeUnits(conf.origin_m),
                unit_converter->ConvertDistanceToLatticeUnits(conf.length_m),
                unit_converter->ConvertDistanceToLatticeUnits(conf.radius_m),
                unit_converter->ConvertDistanceToLatticeUnits(conf.fadelength_m)
        );
    }

    std::shared_ptr<hemelb::net::IteratedAction> SimBuilder::BuildColloidController() const {
        if (config.HasColloidSection()) {
#ifdef HEMELB_BUILD_COLLOIDS
            timings[reporting::Timers::colloidInitialisation].Start();
                    log::Logger::Log<log::Info, log::Singleton>("Loading Colloid config.");
                    std::string colloidConfigPath = simConfig->GetColloidConfigPath();
                    io::xml::Document xml(colloidConfigPath);

                    log::Logger::Log<log::Info, log::Singleton>("Creating Body Forces.");
                    colloids::BodyForces::InitBodyForces(xml);

                    log::Logger::Log<log::Info, log::Singleton>("Creating Boundary Conditions.");
                    colloids::BoundaryConditions::InitBoundaryConditions(domainData.get(), xml);

                    log::Logger::Log<log::Info, log::Singleton>("Initialising Colloids.");
                    auto colloidController =
                            std::make_shared<colloids::ColloidController>(*domainData,
                                                                                  *simulationState,
                                                                                  readGeometryData,
                                                                                  xml,
                                                                                  latticeBoltzmannModel->GetPropertyCache(),
                                                                                  latticeBoltzmannModel->GetLbmParams(),
                                                                                  fileManager->GetColloidPath(),
                                                                                  ioComms,
                                                                                  timings);
                    timings[reporting::Timers::colloidInitialisation].Stop();
#else
            throw Exception() << "Config contains <colloids> tag when built with HEMELB_BUILD_COLLOIDS=OFF";
#endif
        }
        return {};
    }

//    template <typename F>
//    struct transformer {
//        F func;
//
//        transformer(F&& f) : func{f} {}
//
//        template <template <typename...> class ContainerT, typename V, typename... ARGS>
//        //requires std::invocable<F, V const&>
//        auto operator()(ContainerT<V, ARGS...>const & c) {
//            using R = std::invoke_result_t<F, V const&>;
//            ContainerT<R> ans;
//            std::transform(c.begin(), c.end(), std::back_inserter(ans), func);
//            return ans;
//        }
//    };

    std::shared_ptr<extraction::PropertyActor> SimBuilder::BuildPropertyExtraction(
            std::filesystem::path const& xtr_path,
            std::shared_ptr<lb::SimulationState const> simState,
            std::shared_ptr<extraction::IterableDataSource> dataSource,
            reporting::Timers& timings,
            const net::IOCommunicator& ioComms
    ) const {
        // Copy the output file descriptions
        auto po = config.GetPropertyOutputs();
        // Prepend the extraction dir
        for (auto& p: po) {
            p.filename = xtr_path / p.filename;
        }
        // Create the Actor
        return std::make_shared<extraction::PropertyActor>(
                simState,
                po,
                dataSource,
                timings,
                ioComms
        );
    }
    std::shared_ptr<io::Checkpointer> SimBuilder::BuildCheckpointer(
            std::filesystem::path const& cp_path,
            std::shared_ptr<lb::SimulationState const> simState,
            std::shared_ptr<extraction::IterableDataSource> dataSource,
            reporting::Timers& timings,
            const net::IOCommunicator& ioComms
    ) const {
        if (!config.sim_info.checkpoint.has_value())
            return nullptr;

        auto& cp_info = *config.sim_info.checkpoint;

        auto cp_dir_pattern = cp_path / "%d";
        auto dist_pattern = cp_dir_pattern / "distributions.xtr";

        // Create a checkpoint property extractor.
        //
        // This is just a normal one, but fixed to be whole geometry,
        // only distributions, at double precision.
        auto lpo = [&]() {
            using namespace extraction;
            auto file = PropertyOutputFile{
                    .filename = cp_dir_pattern / "distributions.xtr",
                    .frequency = cp_info.period,
                    .geometry = util::make_clone_ptr<WholeGeometrySelector>(),
                    .fields = {
                            OutputField{
                                    .name = "distributions",
                                    .src = source::Distributions{},
                                    .typecode = distribn_t{0},
                                    .noffsets = 0,
                                    .offset = {},
                            }
                    },
                    .ts_mode = single_timestep_files{}
            };
            return std::make_unique<LocalPropertyOutput>(dataSource, file, ioComms);
        }();

        return std::make_shared<io::Checkpointer>(
                cp_info.period, cp_dir_pattern.native(), std::move(simState), std::move(lpo)
        );
    }

    std::shared_ptr<reporting::Reporter> SimBuilder::BuildReporter(
            PathManager const& fileManager,
            std::vector<reporting::Reportable*>const& reps
    ) const {
        auto reporter = std::make_shared<reporting::Reporter>(fileManager.GetReportPath(),
                                                                  fileManager.GetInputFile());
        for (auto r: reps) {
            reporter->AddReportable(r);
        }
        return reporter;
    }
}