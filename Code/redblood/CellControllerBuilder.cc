// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "redblood/CellControllerBuilder.h"

#include "configuration/PathManager.h"
#include "redblood/CellIO.h"
#include "redblood/FaderCell.h"
#include "redblood/MeshIO.h"

namespace hemelb::redblood {

    namespace {
        struct meshio_maker {
            using ptr = std::unique_ptr<MeshIO>;
            ptr operator()(std::monostate) {
                throw Exception() << "Invalid mesh format!";
                return nullptr;
            }

            ptr operator()(configuration::VTKMeshFormat) {
                return std::make_unique<VTKMeshIO>();
            }
            ptr operator()(configuration::KruegerMeshFormat) {
                log::Logger::Log<log::Warning, log::Singleton>("Krueger format meshes are deprecated, move to VTK when you can.");
                return std::make_unique<KruegerMeshIO>();
            }
        };

        bool validateCellEdgeLengths(const CellBase& cell)
        {
            auto edgeLength = cell.GetAverageEdgeLength();

            // Acceptable average edge length to voxel size ratio is [0.7, 1.3].
            // Note that cell vertices location is given in lattice units, therefore
            // no need to normalise again.
            return (edgeLength >= 0.7 && edgeLength <= 1.3);
        }

        //! Rotates a cell to be aligned with the flow and translates it to the start of the flow extension fade length
        void rotateTranslateCellToFlow(CellBase& cell, const Angle theta,
                                       const Angle phi, const FlowExtension & flowExtension,
                                       LatticePosition const & translation,
                                       util::Matrix3D & rotateToFlow, util::Matrix3D & rotation)
        {
            // Rotate cell to align z axis with given position, and then z axis with flow
            // If phi == 0, then cell symmetry axis is aligned with the flow
            using std::cos;
            using std::sin;
            LatticePosition const z(cos(theta) * sin(phi), sin(theta) * sin(phi), cos(phi));
            rotateToFlow = rotationMatrix(LatticePosition(0, 0, 1), flowExtension.normal);
            rotation = rotateToFlow * rotationMatrix(LatticePosition(0, 0, 1), z);
            cell *= rotation;

            // Figure out size of cell alongst cylinder axis
            auto const barycentre = cell.GetBarycentre();
            auto maxExtent = [barycentre, &flowExtension](LatticePosition const pos)
            {
                return std::max(Dot(pos - barycentre, flowExtension.normal), 0e0);
            };
            auto const maxZ =
                    *std::max_element(cell.GetVertices().begin(),
                                      cell.GetVertices().end(),
                                      [&maxExtent](LatticePosition const &a, LatticePosition const& b)
                                      {
                                          return maxExtent(a) < maxExtent(b);
                                      });
            // Place cell as close as possible to 0 of fade length
            cell += flowExtension.origin
                     + flowExtension.normal * (flowExtension.fadeLength - maxExtent(maxZ)) - barycentre
                     + rotateToFlow * translation;

            // fail if any node outside flow extension
            for (auto const &vertex : cell.GetVertices())
            {
                if (not contains(flowExtension, vertex))
                {
                    HEMELB_CAPTURE(flowExtension.normal);
                    HEMELB_CAPTURE(flowExtension.origin);
                    HEMELB_CAPTURE(flowExtension.radius);
                    HEMELB_CAPTURE(flowExtension.length);
                    HEMELB_CAPTURE(vertex);
                    throw Exception() << "BAD INPUT: Cell not contained within flow extension";
                }
            }
        }
    }

    Node2NodeForce CellControllerBuilder::build_node2node_force(configuration::NodeForceConfig const& conf) const {
        // ensure intensity is converted to lattice units
        auto intensity_lat = overload_visit(
                conf.intensity,
                [&](quantity<double, "Nm"> const & q) {
                    return unit_converter->ConvertToLatticeUnits("Nm", q.value());
                },
                [](quantity<double, "lattice"> const& q) {
                    return q.value();
                }
        );

        if (2e0 * conf.cutoffdist > Dimensionless(Traits<>::Stencil::GetRange()))
        {
            log::Logger::Log<log::Warning, log::Singleton>("Input inconsistency: cell-cell and cell-wall interactions larger then stencil size\n"
                                                           "See issue #586.");
            throw Exception() << "Cell-cell interaction longer that stencil size permits";
        }
        return {intensity_lat, conf.cutoffdist, conf.exponent};
    }

    Cell::Moduli CellControllerBuilder::build_cell_moduli(const configuration::CellModuli &conf) const {
        redblood::Cell::Moduli moduli;
        moduli.bending = unit_converter->ConvertToLatticeUnits("Nm", conf.bending_Nm);
        moduli.surface = conf.surface_lat;
        moduli.volume = conf.volume_lat;
        moduli.dilation = conf.dilation_lat;
        moduli.strain = unit_converter->ConvertToLatticeUnits("N/m", conf.strain_Npm);
        return moduli;
    }

    std::unique_ptr<CellBase> CellControllerBuilder::build_cell(configuration::TemplateCellConfig const& tc_conf) const {
        auto meshio = std::visit(meshio_maker{}, tc_conf.format);
        auto mesh_data = meshio->readFile(tc_conf.mesh_path, true);

        auto scale = unit_converter->ConvertDistanceToLatticeUnits(tc_conf.scale_m);
        std::unique_ptr<Cell> cell;
        if (tc_conf.reference_mesh_path.has_value())
        {
            auto io = std::visit(meshio_maker{}, tc_conf.reference_mesh_format);
            auto reference_mesh_data = io->readFile(tc_conf.reference_mesh_path.value(), true);
            if (mesh_data->facets != reference_mesh_data->facets)
                throw Exception() << "Reference mesh facets do not match for cell " << tc_conf.name;
            if (!(volume(mesh_data->vertices, reference_mesh_data->facets) > 0.0))
                throw Exception() << "Reference mesh volume calculation not positive for cell " << tc_conf.name;

            cell = std::make_unique<Cell>(mesh_data->vertices, reference_mesh_data, scale, tc_conf.name);
        } else {
            cell = std::make_unique<Cell>(mesh_data->vertices, Mesh(mesh_data), scale, tc_conf.name);
        }
        *cell *= scale;
        cell->moduli = build_cell_moduli(tc_conf.moduli);

        return cell;
    }

    std::shared_ptr<TemplateCellContainer> CellControllerBuilder::build_template_cells(
            const std::map<std::string, configuration::TemplateCellConfig> &conf,
            CountedIoletView const& inlets,
            CountedIoletView const&  outlets) const {
        auto result = std::make_shared<TemplateCellContainer>();

        // First, concatenate all flow extensions, if they exist.
        // Needs to be a shared pointer so FaderCell can share ownership.
        auto flowExtensions = [&]() {
            using VFE = std::vector<FlowExtension>;
            VFE ans;
            auto add_fext = [&ans] (CountedIoletView const& iolets) {
                for (int i = 0; i < iolets.GetCount(); ++i) {
                    auto iolet = iolets.GetIolet(i);
                    if (auto fext = iolet->GetFlowExtension())
                        ans.push_back(*fext);
                }
            };
            add_fext(inlets);
            add_fext(outlets);

            if (ans.empty()) {
                return std::shared_ptr<VFE>(nullptr);
            } else {
                return std::make_shared<VFE>(std::move(ans));
            }
        }();

        // Construct the templates
        for (auto [name, tc_conf]: conf) {
            auto cell = build_cell(tc_conf);
            if (!validateCellEdgeLengths(*cell))
            {
                log::Logger::Log<log::Critical, log::Singleton>("Average edge length in cell mesh not in [0.7, 1.3]");
            }
            if (flowExtensions)
            {
                auto fader = FaderCell(std::move(cell), flowExtensions).clone();
                cell = std::move(fader);
            }
            result->emplace(name, std::shared_ptr<CellBase>(cell.release()));
        }
        return result;
    }

    CompositeRBCInserter CellControllerBuilder::build_single_inlet_rbc_inserter(
            std::vector<configuration::CellInserterConfig> const& ci_confs,
            lb::InOutLet const& inlet,
            TemplateCellContainer const& templateCells
    ) const {
        CompositeRBCInserter composite;
        for (auto& conf: ci_confs) {
            // Clone, since we rotate and shift below.
            auto cell = templateCells.at(conf.template_name)->clone();

            auto flowExtension = inlet.GetFlowExtension();

            // Rotate cell to align z axis with given position, and then z axis with flow
            // If phi == 0, then cell symmetry axis is aligned with the flow
            auto const theta = conf.theta_rad;
            auto const phi = conf.phi_rad;
            LatticeDisplacement  trans = unit_converter->ConvertDisplacementToLatticeUnits(conf.translation_m);

            util::Matrix3D rotateToFlow, rotation;
            rotateTranslateCellToFlow(*cell,
                                      theta,
                                      phi,
                                      *flowExtension,
                                      trans,
                                      rotateToFlow,
                                      rotation);

            // Drops first cell when time reaches offset, and then every deltaTime thereafter.
            // Note: c++14 will allow more complex captures. Until then, we will need to create
            // semi-local lambda variables on the stack as shared pointers. Where semi-local means the
            // variables should live as long as the lambda. But longer than a single call.
            auto const offset = overload_visit(
                    conf.offset,
                    [&](quantity<double, "s"> const & q){
                        return unit_converter->ConvertTimeToLatticeUnits(q.value());
                    },
                    [](quantity<double, "lattice"> const& q) {
                        return q.value();
                    }
            );
            auto const timeStep = unit_converter->ConvertTimeToLatticeUnits(conf.drop_period_s);
            auto const dt = unit_converter->ConvertTimeToLatticeUnits(conf.dt_s);
            auto time = std::make_shared<LatticeTime>(timeStep - 1e0
                                                      + std::numeric_limits<LatticeTime>::epsilon() - offset);

            std::minstd_rand randomGenerator(conf.seed);
            std::uniform_real_distribution<double> uniformDistribution(-1.0, 1.0);

            auto condition = [time, timeStep, dt, uniformDistribution, randomGenerator]() mutable
            {
                *time += 1e0;
                if(*time >= timeStep)
                {
                    *time -= timeStep + dt * uniformDistribution(randomGenerator);
                    return true;
                }
                return false;
            };

            auto const dx = unit_converter->ConvertDistanceToLatticeUnits(conf.dx_m);
            auto const dy = unit_converter->ConvertDistanceToLatticeUnits(conf.dy_m);

            composite.AddInserter(std::static_pointer_cast<RBCInserter>(std::make_shared<
                    RBCInserterWithPerturbation>(condition,
                                                 std::move(cell),
                                                 rotation,
                                                 conf.dtheta_rad,
                                                 conf.dphi_rad,
                                                 rotateToFlow * LatticePosition(dx, 0, 0),
                                                 rotateToFlow * LatticePosition(0, dy, 0),
                                                 conf.seed)));
        }

        return composite;
    }

    std::function<void(CellInserter const&)>
    CellControllerBuilder::build_cell_inserters(const std::vector<configuration::IoletConfig> &inlet_confs,
                                                CountedIoletView const& inlets,
                                                const TemplateCellContainer &templateCells) const {
        // Check if we have any cell insertion action going on at all
        //auto inlet = findFirstInletWithCellInsertion(inletsNode.GetChildOrThrow("inlet"));
        std::vector<std::function<void(CellInserter const&)>> results;
        for (auto [i, inlet_conf]: util::enumerate(inlet_confs)) {
            auto ci_confs =  std::visit([](auto&& ic) -> std::vector<configuration::CellInserterConfig> {
                using T = std::decay_t<decltype(ic)>;
                if constexpr (std::is_same_v<T, std::monostate>) {
                    return {};
                } else {
                    if (ic.flow_extension.has_value())
                        return ic.cell_inserters;
                    return {};
                }
            }, inlet_conf);
            if (!ci_confs.empty())
                results.emplace_back(build_single_inlet_rbc_inserter(ci_confs, *inlets.GetIolet(i), templateCells));
        }

        // do all insertion functions in one go
        if (results.size() > 1)
        {
            // go to shared pointer to avoid copies
            std::shared_ptr<decltype(results)> functions(new decltype(results)(std::move(results)));
            // return a lambda that loops over functions
            return [functions](CellInserter const& inserter)
            {
                for(auto const & func: *functions)
                {
                    func(inserter);
                }
            };
        }
        // return null if no insertion, and just the function if only one
        return results.size() == 0 ?
               nullptr :
               results.front();
    }

    std::vector<FlowExtension> CellControllerBuilder::build_outlets(
            std::vector<configuration::IoletConfig> const& inlet_confs,
            CountedIoletView const& inlets,
            CountedIoletView const& outlets
    ) const {
        using VFE = std::vector<FlowExtension>;
        VFE result;

        // Transforms them to cell outlets: should start somewhere near the end of fadelength
        auto trans = [] (FlowExtension flowExt) {
            auto length = flowExt.length - flowExt.fadeLength;
            flowExt.origin += flowExt.normal * (flowExt.length - length);
            flowExt.length = length;
            flowExt.fadeLength = length;
            return flowExt;
        };

        for (int i = 0; i < outlets.GetCount(); ++i) {
            auto& outlet = *outlets.GetIolet(i);
            if (auto fext = outlet.GetFlowExtension())
                result.push_back(trans(*fext));
        }

        for (int i = 0; i < inlets.GetCount(); ++i) {
            auto& inlet = *inlets.GetIolet(i);
            bool has_flow_ext_but_no_inserters = std::visit([](auto&& icv) {
                using T = std::decay_t<decltype(icv)>;
                if constexpr (std::is_same_v<T, std::monostate>) {
                    return false;
                } else {
                    return icv.flow_extension.has_value() && icv.cell_inserters.empty();
                }
            }, inlet_confs[i]);
            if (has_flow_ext_but_no_inserters)
                result.push_back(trans(*inlet.GetFlowExtension()));
        }

        return result;
    }

    CellChangeListener CellControllerBuilder::build_full_cell_output(
            configuration::CellOutputConfig const& conf,
            std::shared_ptr<lb::SimulationState const> simState,
            std::shared_ptr<configuration::PathManager const> fileManager,
            net::IOCommunicator const& ioComms
    ) const {
        auto conv = conf.physical_units ? unit_converter : nullptr;
        return CellVtkOutput{conf.output_period, conv, simState, fileManager, ioComms};
    }

    CellChangeListener CellControllerBuilder::build_summary_cell_output(
            configuration::CellOutputConfig const& conf,
            std::shared_ptr<lb::SimulationState const> simState,
            std::shared_ptr<configuration::PathManager const> fileManager,
            net::IOCommunicator const& ioComms
    ) const {
        auto conv = conf.physical_units ? unit_converter : nullptr;
        return CellBarycentreOutput{conf.output_period, conv, simState, fileManager, ioComms};
    }

}