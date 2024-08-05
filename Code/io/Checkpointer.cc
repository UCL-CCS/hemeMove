// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "io/Checkpointer.h"

#include "configuration/SimConfig.h"
#include "configuration/SimConfigWriter.h"
#include "extraction/LocalPropertyOutput.h"

namespace hemelb::io {
    using namespace configuration;

    Checkpointer::Checkpointer(
            LatticeTimeStep p,
            std::string_view pattern,
            std::shared_ptr<lb::SimulationState const> ss,
            std::unique_ptr<extraction::LocalPropertyOutput> dist_writer
    ) :
            period(p),
            out_dir_pattern(pattern),
            simulationState(std::move(ss)),
            distribution_writer(std::move(dist_writer))
    {}

    Checkpointer::~Checkpointer() = default;

    bool Checkpointer::ShouldWrite() const {
        return (simulationState->GetTimeStep() % period) == 0;
    }

    void Checkpointer::Write(SimConfig conf) const {
        auto const dir = path(out_dir_pattern.Format(simulationState->GetTimeStep(), simulationState->GetEndTimeStep()));
        std::filesystem::create_directory(dir);
        auto [dist_path, off_path] = WriteDists();

        auto current_time = simulationState->GetTimeStep();
        conf.initial_condition = CheckpointIC(current_time, dist_path, std::make_optional(std::move(off_path)));

        // Update RBC conf with seeds etc
        if (distribution_writer->OnIORank()) {
            std::filesystem::path const xml_path = dir / "restart.xml";
            WriteConfig(conf, xml_path);
        }
    }

    auto Checkpointer::WriteDists() const -> std::pair<path, path> {
        auto dist_path = distribution_writer->Write(simulationState->GetTimeStep(), simulationState->GetEndTimeStep()).value();
        auto off_path = distribution_writer->GetOffsetFileName();
        return {dist_path, off_path};
    }

    void Checkpointer::WriteConfig(SimConfig const& conf, path const& out_xml) const {
        auto w = SimConfigWriter(out_xml);
        w.Write(conf);
    }
}