// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "configuration/SimConfig.h"
#include "configuration/SimConfigReader.h"

namespace hemelb::configuration
{

    // Base IC
    ICConfigBase::ICConfigBase(std::optional<LatticeTimeStep> t) : t0(t) {
    }

    // Uniform equilibrium IC
    EquilibriumIC::EquilibriumIC(std::optional<LatticeTimeStep> t, PhysicalPressure p) : ICConfigBase(t), p_Pa(p), v_ms(0.0) {
    }

    EquilibriumIC::EquilibriumIC(std::optional<LatticeTimeStep> t, PhysicalPressure p, const PhysicalVelocity& v) : ICConfigBase(t), p_Pa(p), v_ms(v) {
    }

    // checkpoint IC
    CheckpointIC::CheckpointIC(std::optional<LatticeTimeStep> t, std::filesystem::path cp, std::optional<std::filesystem::path> maybeOff)
            : ICConfigBase(t), cpFile(std::move(cp)), maybeOffFile(std::move(maybeOff)) {
    }

    SimConfig SimConfig::New(const path& path)
    {
        auto reader = SimConfigReader(path);
        return reader.Read();
    }

    const MonitoringConfig& SimConfig::GetMonitoringConfiguration() const
    {
      return monitoringConfig;
    }

}
