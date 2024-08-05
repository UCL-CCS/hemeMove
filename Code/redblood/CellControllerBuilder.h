// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_REDBLOOD_CELLCONTROLLERBUILDER_H
#define HEMELB_REDBLOOD_CELLCONTROLLERBUILDER_H

#include <memory>
#include <optional>

#include "configuration/SimConfig.h"
#include "net/IteratedAction.h"
#include "util/UnitConverter.h"

#include "redblood/types.h"
#include "redblood/CellController.h"
#include "redblood/RBCInserter.h"

namespace hemelb::configuration { class PathManager; }

namespace hemelb::redblood {

    //! For testing purposes, we don't want to have to create the
    //! BoundaryValues objects to get at iolets, so create a simple ,
    //! type-erased wrapper which accesses them with the same style.
    class CountedIoletView {
        using Getter = std::function<lb::InOutLet const* (unsigned)>;
        using Counter = std::function<unsigned()>;

        Counter counter;
        Getter getter;

    public:
        // Counter must give the number of items and Getter gets one by index.
        template <typename CounterT, typename GetterT>
        CountedIoletView(CounterT && c, GetterT&& g) : counter{c}, getter{g} {

        }
        // Invoke the counter
        inline unsigned GetCount() const {
            return counter();
        }
        // Invoke the getter
        inline lb::InOutLet const* GetIolet(unsigned i) const {
            return getter(i);
        }
    };


    // This class performs much the same role as SimBuilder, but for the RBC code only.
    // The main entry point is the `build` member function template below.
    // TODO: document requirements on what of the main simulation must have been created first.
    class CellControllerBuilder {
        std::shared_ptr<util::UnitConverter const> unit_converter;
    public:
        inline explicit CellControllerBuilder(std::shared_ptr<util::UnitConverter const> uc)
                : unit_converter(std::move(uc)) {
        }

        std::shared_ptr<TemplateCellContainer> build_template_cells(
                std::map<std::string, configuration::TemplateCellConfig> const& conf,
                CountedIoletView const& inlets,
                CountedIoletView const& outlets
        ) const;
        Cell::Moduli build_cell_moduli(configuration::CellModuli const& conf) const;
        std::unique_ptr<CellBase> build_cell(configuration::TemplateCellConfig const& tc_conf) const;
        Node2NodeForce build_node2node_force(configuration::NodeForceConfig const&) const;
        CompositeRBCInserter build_single_inlet_rbc_inserter(
                std::vector<configuration::CellInserterConfig> const& ci_confs,
                lb::InOutLet const& inlet,
                TemplateCellContainer const& templateCells
        ) const;

        std::function<void(CellInserter const&)> build_cell_inserters(
                std::vector<configuration::IoletConfig> const& inlet_confs,
                CountedIoletView const& inlets,
                TemplateCellContainer const& templateCells
        ) const;

        std::vector<FlowExtension> build_outlets(
                std::vector<configuration::IoletConfig> const& inlet_confs,
                CountedIoletView const& inlets,
                CountedIoletView const& outlets
        ) const;
        CellChangeListener build_full_cell_output(
                configuration::CellOutputConfig const& conf,
                std::shared_ptr<lb::SimulationState const> simState,
                std::shared_ptr<configuration::PathManager const> fileManager,
                net::IOCommunicator const& ioComms
        ) const;
        CellChangeListener build_summary_cell_output(
                configuration::CellOutputConfig const& conf,
                std::shared_ptr<lb::SimulationState const> simState,
                std::shared_ptr<configuration::PathManager const> fileManager,
                net::IOCommunicator const& ioComms
        ) const;

        template <typename Traits>
        std::shared_ptr<net::IteratedAction> build(configuration::SimConfig const& config,
                                                   reporting::Timers& timings,
                                                   net::IOCommunicator const& ioComms,
                                                   std::shared_ptr<geometry::FieldData> fieldData,
                                                   CountedIoletView const& inlets,
                                                   CountedIoletView const& outlets,
                                                   std::shared_ptr<lb::SimulationState const> simState,
                                                   std::shared_ptr<configuration::PathManager const> fileManager
        ) {
            log::Logger::Log<log::Info, log::Singleton>("Initialising RBCs.");
            timings.cellInitialisation().Start();
            auto& rbcConfig = config.GetRBCConfig();
            CellContainer cells;
            using Controller = CellController<Traits>;
            auto meshes = build_template_cells(rbcConfig.meshes, inlets, outlets);
            auto controller = std::make_shared<Controller>(
                    *fieldData,
                    cells,
                    meshes,
                    timings,
                    rbcConfig.boxSize,
                    build_node2node_force(rbcConfig.cell2cell),
                    build_node2node_force(rbcConfig.cell2wall),
                    ioComms);

            controller->SetCellInsertion(build_cell_inserters(config.GetInlets(), inlets, *meshes));

            controller->SetOutlets(build_outlets(config.GetInlets(), inlets, outlets));

            // Lambda to DRY
            auto maybe_add_cell_output = [&](std::optional<configuration::CellOutputConfig> const& maybe_conf, auto factory) {
                if (maybe_conf) {
                    controller->AddCellChangeListener(factory(
                            maybe_conf.value(),
                            simState,
                            fileManager,
                            ioComms
                    ));
                }
            };
            maybe_add_cell_output(rbcConfig.full_output,
                                  [this](auto... args){ return this->build_full_cell_output(args...); });
            maybe_add_cell_output(rbcConfig.summary_output,
                                  [this](auto... args){ return this->build_summary_cell_output(args...); });

            timings.cellInitialisation().Stop();
            return controller;
        }
    };
}
#endif