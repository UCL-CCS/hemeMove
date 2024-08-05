// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "configuration/SimConfigWriter.h"

#include <optional>
#include <sstream>
#include <type_traits>
#include <utility>
#include <variant>

#include "constants.h"
#include "Exception.h"
#include "hassert.h"
#include "configuration/MonitoringConfig.h"
#include "extraction/GeometrySelectors.h"
#include "io/xml.h"
#include "util/variant.h"

namespace hemelb::configuration {
    using namespace io::xml;

    namespace {
        template<typename T>
        void SetDimensionalValue(Element el, char const *unit, T const &value) {
            el.SetAttribute("units", unit);
            el.SetAttribute("value", value);
        }

        template<typename T>
        void AddChildDimensionalValue(Element parent, char const *name, char const *unit, T const &value) {
            auto el = parent.AddChild(name);
            SetDimensionalValue(el, unit, value);
        }

        template<typename CStr, typename... CStrs>
        std::optional<Element> GetElemMaybe(Element el, CStr nextPath, CStrs... restPath) {
            auto nextEl = el.GetChildOrNull(nextPath);
            if (nextEl == Element::Missing())
                return std::nullopt;
            if constexpr (sizeof...(CStrs) == 0) {
                return nextEl;
            } else {
                return GetElemMaybe(nextEl, restPath...);
            }
        }

        template<typename... CStr>
        std::optional<Element> GetElemMaybe(std::shared_ptr<const Document> const &doc, CStr... path) {
            return GetElemMaybe(doc->GetRoot(), path...);
        }
    }

//    template <typename... CStrs>
//    bool SimConfigWriter::InOriginal(CStrs... path) const {
//        return GetElemMaybe(originalXml->GetRoot(), path...).has_value();
//    }

    namespace {
        // Make a Document that ensures floats are written with full precision.
        auto MakeXmlDoc() {
            return std::make_unique<Document>([]() {
                std::ostringstream s;
                s << std::hexfloat;
                return s;
            });
        }
    }

    SimConfigWriter::SimConfigWriter(path p)
            : outputXmlPath(std::move(p)), outputXml(MakeXmlDoc())
    {
    }

    void SimConfigWriter::AddPathChild(Element& parent,
                                       const path& p) const {
        auto path = parent.AddChild("path");
        path.SetAttribute("value", FullPathToRelPath(p));
    }

    std::string SimConfigWriter::FullPathToRelPath(const SimConfigWriter::path& p) const {
        auto out_dir = outputXmlPath.parent_path();
        return p.lexically_relative(out_dir);
    }


    void SimConfigWriter::Write(SimConfig const& conf) {
        auto root = outputXml->AddChild("hemelbsettings");
        root.SetAttribute("version", 6U);
        DoIOForSimulation(conf.GetSimInfo());
        DoIOForGeometry(conf.GetDataFilePath());
        DoIOForInitialConditions(conf.GetInitialCondition());
        DoIOForInOutlets("inlet", conf.GetInlets());
        DoIOForInOutlets("outlet", conf.GetOutlets());
        DoIOForProperties(conf.GetPropertyOutputs());
        DoIOForMonitoring(conf.GetMonitoringConfiguration());
        // DoIOForRBCs

        outputXml->SaveFile(outputXmlPath);
        // Reset to a new empty doc.
        outputXml = MakeXmlDoc();
    }

    void SimConfigWriter::DoIOForSimulation(const GlobalSimInfo& sim_info) {
        Element outSimEl = outputXml->GetRoot().AddChild("simulation");

        AddChildDimensionalValue(outSimEl, "step_length", "s", sim_info.time.step_s);

        // Required element
        // <steps value="unsigned" units="lattice />
        AddChildDimensionalValue(outSimEl, "steps", "lattice", sim_info.time.total_steps);

        // Optional element
        // <extra_warmup_steps value="unsigned" units="lattice" />
        if (auto warmup = sim_info.time.warmup_steps; warmup > 0)
            AddChildDimensionalValue(outSimEl, "extra_warmup_steps", "lattice", warmup);

        // Required element
        // <voxel_size value="float" units="m" />
        AddChildDimensionalValue(outSimEl, "voxel_size", "m", sim_info.space.step_m);
        // Required element
        // <origin value="(x,y,z)" units="m" />
        AddChildDimensionalValue(outSimEl, "origin", "m", sim_info.space.geometry_origin_m);

        // Optional element
        // <fluid_density value="float" units="kg/m3" />
        if (auto density = sim_info.fluid.density_kgm3; density != DEFAULT_FLUID_DENSITY_Kg_per_m3)
            AddChildDimensionalValue(outSimEl, "fluid_density", "kg/m3", density);
        // Optional element
        // <fluid_viscosity value="float" units="Pa.s" />
        if (auto visc = sim_info.fluid.viscosity_Pas; visc != DEFAULT_FLUID_VISCOSITY_Pas)
            AddChildDimensionalValue(outSimEl, "fluid_viscosity", "Pa.s", visc);

        // Optional element
        // <reference_pressure value="float" units="Pa" />
        if (auto ref_p = sim_info.fluid.reference_pressure_Pa; ref_p != 0.0)
            AddChildDimensionalValue(outSimEl, "reference_pressure", "Pa", ref_p);


        if (sim_info.checkpoint.has_value()) {
            auto cpEl = outSimEl.AddChild("checkpoint");
            cpEl.SetAttribute("period", sim_info.checkpoint->period);
        }
    }

    void SimConfigWriter::DoIOForGeometry(const path& p) {
        auto rel_path = FullPathToRelPath(p);
        auto df_el = outputXml->GetRoot().AddChild("geometry").AddChild("datafile");
        df_el.SetAttribute("path", rel_path);
    }

    void SimConfigWriter::DoIOForBaseInOutlet(Element outEl, IoletConfigBase const& ioletConf) const {
        AddChildDimensionalValue(outEl, "position", "m", ioletConf.position);
        AddChildDimensionalValue(outEl, "normal", "dimensionless", ioletConf.normal);

        if (ioletConf.flow_extension.has_value()) {
            auto& fe = ioletConf.flow_extension.value();
            auto flowEl = outEl.AddChild("flowextension");
            AddChildDimensionalValue(flowEl, "length", "m", fe.length_m);
            AddChildDimensionalValue(flowEl, "radius", "m", fe.radius_m);
        }
    }

    void SimConfigWriter::DoIOForInOutlets(std::string type, const std::vector<IoletConfig>& iolets) {
        auto plural = type + "s";

        auto dst_plural = outputXml->GetRoot().AddChild(plural.c_str());

        for (auto& iolet_conf: iolets) {
            auto dst_iolet_el = dst_plural.AddChild(type.c_str());
            std::visit([&](auto const& _) {
                if constexpr(std::is_same_v<std::decay_t<decltype(_)>, std::monostate>) {
                    throw (Exception() << "invalid iolet config");
                } else {
                    DoIOForBaseInOutlet(dst_iolet_el, _);
                }
            }, iolet_conf);
            overload_visit(iolet_conf,
                           [](std::monostate const &) {
                               throw (Exception() << "Invalid iolet config");
                           },
                           [&](CosinePressureIoletConfig const& _) {
                               DoIOForCosinePressureInOutlet(dst_iolet_el, _);
                           },
                           [&](FilePressureIoletConfig const& _) {
                               DoIOForFilePressureInOutlet(dst_iolet_el, _);
                           },
                           [&](MultiscalePressureIoletConfig const& _) {
                               DoIOForMultiscalePressureInOutlet(dst_iolet_el, _);
                           },
                           [&](ParabolicVelocityIoletConfig const& _) {
                               DoIOForParabolicVelocityInOutlet(dst_iolet_el, _);
                           },
                           [&](WomersleyVelocityIoletConfig const& _) {
                               DoIOForWomersleyVelocityInOutlet(dst_iolet_el, _);
                           },
                           [&](FileVelocityIoletConfig const& _) {
                               DoIOForFileVelocityInOutlet(dst_iolet_el, _);
                           }
            );
        }
    }
    Element MakeCondition(Element& iolet, char const* type, char const* subtype) {
        auto condition = iolet.AddChild("condition");
        condition.SetAttribute("type", type);
        condition.SetAttribute("subtype", subtype);
        return condition;
    }

    void SimConfigWriter::DoIOForCosinePressureInOutlet(SimConfigWriter::Element& dest,
                                                        CosinePressureIoletConfig const& conf) const {
        auto condition = MakeCondition(dest, "pressure", "cosine");
        AddChildDimensionalValue(condition, "amplitude", "Pa", conf.amp_Pa);
        AddChildDimensionalValue(condition, "mean", "Pa", conf.mean_Pa);
        AddChildDimensionalValue(condition, "phase", "rad", conf.phase_rad);
        AddChildDimensionalValue(condition, "period", "s", conf.period_s);
    }

    void SimConfigWriter::DoIOForFilePressureInOutlet(SimConfigWriter::Element& dest,
                                                      FilePressureIoletConfig const& conf) const {
        auto condition = MakeCondition(dest, "pressure", "file");
        AddPathChild(condition, conf.file_path);
    }

    void SimConfigWriter::DoIOForMultiscalePressureInOutlet(SimConfigWriter::Element& dest,
                                                            MultiscalePressureIoletConfig const& conf) const {
        auto condition = MakeCondition(dest, "pressure", "multiscale");
        AddChildDimensionalValue(condition, "pressure", "Pa", conf.pressure_reference_Pa);
        AddChildDimensionalValue(dest, "velocity", "m/s", conf.velocity_reference_ms);
        auto label = condition.AddChild("label");
        label.SetAttribute("value", conf.label);
    }

    void SimConfigWriter::DoIOForParabolicVelocityInOutlet(SimConfigWriter::Element& dest,
                                                           ParabolicVelocityIoletConfig const& conf) const {
        auto condition = MakeCondition(dest, "velocity", "parabolic");
        AddChildDimensionalValue(dest, "radius", "m", conf.radius_m);
        AddChildDimensionalValue(dest, "maximum", "m/s", conf.max_speed_ms);
    }

    void SimConfigWriter::DoIOForWomersleyVelocityInOutlet(SimConfigWriter::Element& dest,
                                                           WomersleyVelocityIoletConfig const& conf) const {
        auto condition = MakeCondition(dest, "velocity", "womersley");
        AddChildDimensionalValue(dest, "radius", "m", conf.radius_m);
        AddChildDimensionalValue(dest, "pressure_gradient_amplitude", "Pa/m", conf.pgrad_amp_Pam);
        AddChildDimensionalValue(dest, "period", "s", conf.period_s);
        AddChildDimensionalValue(dest, "womersley_number", "dimensionless", conf.womersley);
    }

    void SimConfigWriter::DoIOForFileVelocityInOutlet(SimConfigWriter::Element& dest,
                                                      FileVelocityIoletConfig const& conf) const {
        auto condition = MakeCondition(dest, "velocity", "file");
        AddPathChild(condition, conf.file_path);
    }

    void SimConfigWriter::DoIOForProperties(std::vector<extraction::PropertyOutputFile> const& outputs) const {
        using namespace extraction;

        if (outputs.empty())
            return;
        auto prop_el = outputXml->GetRoot().AddChild("properties");
        for (auto& po: outputs) {
            auto po_el = prop_el.AddChild("propertyoutput");

            auto mode = overload_visit(po.ts_mode,
                                       [](multi_timestep_file) {
                                           return "multi";
                                       },
                                       [](single_timestep_files) {
                                           return "single";
                                       }
            );
            po_el.SetAttribute("timestep_mode", mode);
            po_el.SetAttribute("file", po.filename.c_str());
            po_el.SetAttribute("period", po.frequency);

            GeometrySelector* gmy = po.geometry.get();
            auto gmy_el = po_el.AddChild("geometry");
            if (auto whole = dynamic_cast<WholeGeometrySelector*>(gmy)) {
                gmy_el.SetAttribute("type", "whole");
            } else if (auto surface = dynamic_cast<SurfacePointSelector*>(gmy)) {
                gmy_el.SetAttribute("type", "surface");
            } else if (auto plane = dynamic_cast<PlaneGeometrySelector*>(gmy)) {
                gmy_el.SetAttribute("type", "plane");
                AddChildDimensionalValue(gmy_el, "point", "m", plane->GetPoint());
                AddChildDimensionalValue(gmy_el, "normal", "dimensionless", plane->GetNormal());
                if (auto r = plane->GetRadius(); r > 0.0f) {
                    AddChildDimensionalValue(gmy_el, "radius", "m", r);
                }
            } else if (auto line = dynamic_cast<StraightLineGeometrySelector*>(gmy)) {
                gmy_el.SetAttribute("type", "line");
                AddChildDimensionalValue(gmy_el, "point", "m", line->GetEndpoint1());
                AddChildDimensionalValue(gmy_el, "point", "m", line->GetEndpoint2());
            } else if (auto surf_pt = dynamic_cast<SurfacePointSelector*>(gmy)) {
                gmy_el.SetAttribute("type", "surfacepoint");
                AddChildDimensionalValue(gmy_el, "point", "m", surf_pt->GetPoint());
            } else {
                throw (Exception() << "unknown type for property extraction geometry");
            }

            for (auto& f: po.fields) {
                auto field_el = po_el.AddChild("field");
                field_el.SetAttribute(
                        "type",
                        overload_visit(f.src,
                                       [](source::Pressure) {
                                           return "pressure";
                                       },
                                       [](source::Velocity) {
                                           return "velocity";
                                       },
                                       [](source::VonMisesStress) {
                                           return "vonmisesstress";
                                       },
                                       [](source::ShearStress) {
                                           return "shearstress";
                                       },
                                       [](source::ShearRate) {
                                           return "shearrate";
                                       },
                                       [](source::StressTensor) {
                                           return "stresstensor";
                                       },
                                       [](source::Traction) {
                                           return "traction";
                                       },
                                       [](source::TangentialProjectionTraction) {
                                           return "tangentialprojectiontraction";
                                       },
                                       [](source::Distributions) {
                                           return "distributions";
                                       },
                                       [](source::MpiRank) {
                                           return "mpirank";
                                       }
                        )
                );
                field_el.SetAttribute("name", f.name);
            }
        }
    }

    void SimConfigWriter::DoIOForInitialConditions(ICConfig const& ic_conf) {
        auto ic_el = outputXml->GetRoot().AddChild("initialconditions");

        auto set_time_maybe = [&](ICConfigBase const& base_conf) {
            if (base_conf.t0.has_value())
                AddChildDimensionalValue(ic_el, "time", "lattice", base_conf.t0.value());
        };

        overload_visit(ic_conf,
                       [](std::monostate const &) {
                           throw (Exception() << "Invalid initial condition config");
                       },
                       [&](EquilibriumIC const& _) {
                           set_time_maybe(_);
                           auto p_el = ic_el.AddChild("pressure");
                           AddChildDimensionalValue(p_el, "uniform", "Pa", _.p_Pa);
                       },
                       [&](CheckpointIC const& _) {
                           set_time_maybe(_);
                           auto cp_el = ic_el.AddChild("checkpoint");
                           cp_el.SetAttribute("file", FullPathToRelPath(_.cpFile));
                           if (_.maybeOffFile.has_value())
                               cp_el.SetAttribute("offsets", FullPathToRelPath(_.maybeOffFile.value()));
                       }
        );
    }

    void SimConfigWriter::DoIOForMonitoring(const hemelb::configuration::MonitoringConfig &mon_conf) {
        auto mon_el = Element::Missing();
        auto ensure_monitoring_el = [&] () {
            if (mon_el == Element::Missing())
                mon_el = outputXml->GetRoot().AddChild("monitoring");
            return mon_el;
        };

        if (mon_conf.doConvergenceCheck) {
            auto el = ensure_monitoring_el();
            auto conv_el = el.AddChild("steady_flow_convergence");
            conv_el.SetAttribute("tolerance", mon_conf.convergenceRelativeTolerance);
            conv_el.SetAttribute("terminate", mon_conf.convergenceTerminate ? "true" : "false");

            HASSERT(std::holds_alternative<extraction::source::Velocity>(mon_conf.convergenceVariable));
            auto crit_el = conv_el.AddChild("criterion");
            crit_el.SetAttribute("type", "velocity");
            SetDimensionalValue(crit_el, "m/s", mon_conf.convergenceReferenceValue);
        }

        if (mon_conf.doIncompressibilityCheck) {
            auto el = ensure_monitoring_el();
            auto incomp_el = el.AddChild("incompressibility");
        }
    }
}