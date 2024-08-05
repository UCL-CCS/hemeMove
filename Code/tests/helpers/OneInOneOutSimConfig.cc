// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "tests/helpers/OneInOneOutSimConfig.h"

namespace hemelb::tests::helpers
{
    using namespace configuration;
    using namespace io::xml;

    SimConfig OneInOneOutSimConfigReader::Read() const {
        Document doc;
        auto root = doc.AddChild("hemelbsettings");
        root.SetAttribute("version", 6);

        root.AddChild("simulation");
        root.AddChild("geometry").AddChild("datafile").SetAttribute("path", "");
        root.AddChild("initialconditions");
        root.AddChild("inlets");
        root.AddChild("outlets");

        return DoIO(root);
    }

    GlobalSimInfo OneInOneOutSimConfigReader::DoIOForSimulation(const io::xml::Element simEl) const {
        GlobalSimInfo sim_info;
        sim_info.time = {10000, 0, 60.0 / (70.0 * 1000.0)};
        sim_info.space = {0.01, util::Vector3D<PhysicalDistance>::Zero()};
        sim_info.fluid = {DEFAULT_FLUID_DENSITY_Kg_per_m3, DEFAULT_FLUID_VISCOSITY_Pas, 80.0};
        return sim_info;
    }

    ICConfig OneInOneOutSimConfigReader::DoIOForInitialConditions(io::xml::Element parent) const {
        return {};
    }

    std::vector<IoletConfig>
    OneInOneOutSimConfigReader::DoIOForInOutlets(const GlobalSimInfo &sim_info,
                                                 const io::xml::Element xmlNode) const {
        if (xmlNode.GetName() == "inlets") {
            configuration::CosinePressureIoletConfig inlet;
            inlet.amp_Pa = 1.0;
            inlet.mean_Pa = 80.0;
            inlet.phase_rad = PI;
            inlet.period_s = 60.0 / 70.0;
            inlet.normal = util::Vector3D<Dimensionless>(-3, 4, -9);
            return {inlet};
        }
        if (xmlNode.GetName() == "outlets") {
            configuration::CosinePressureIoletConfig outlet;
            outlet.amp_Pa = 0.0;
            outlet.mean_Pa = 80.0;
            outlet.phase_rad = 0.0;
            outlet.period_s = 60.0 / 70.0;
            outlet.normal = util::Vector3D<Dimensionless>(2, -1, 4);
            return {outlet};
        }
        throw (Exception() << "Not inlets or outlets?");
    }

    OneInOneOutSimConfigReader::OneInOneOutSimConfigReader()
            : configuration::SimConfigReader(std::filesystem::path())
    {

    }

    void OneInOneOutSimConfigReader::CheckIoletMatchesCMake(const io::xml::Element& ioletEl,
                                                            std::string_view requiredBC) const
    {
    }

}
