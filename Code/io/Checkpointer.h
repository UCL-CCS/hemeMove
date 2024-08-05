// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_IO_CHECKPOINTER_H
#define HEMELB_IO_CHECKPOINTER_H

#include <filesystem>
#include <memory>

#include "configuration/SimConfig.h"
#include "io/TimePattern.h"
#include "lb/SimulationState.h"
#include "net/IteratedAction.h"

namespace hemelb::extraction { class LocalPropertyOutput; }

namespace hemelb::io {

    // Control the writing of a complete representation of a simulation state in order to restart.
    class Checkpointer {
        using path = std::filesystem::path;

        LatticeTimeStep period;
        // Path to each individual checkpoint's pattern, with a single "%d" for the time
        io::TimePattern out_dir_pattern;
        std::shared_ptr<lb::SimulationState const> simulationState;
        // Delegate writing of the distribution array to property extraction
        std::unique_ptr<extraction::LocalPropertyOutput> distribution_writer;

    public:
        Checkpointer(
                LatticeTimeStep period,
                std::string_view cp_file_pattern,
                std::shared_ptr<lb::SimulationState const> ss,
                std::unique_ptr<extraction::LocalPropertyOutput> dist_writer
        );

        // Even though trivial, we declare to avoid requiring the definition of LocalPropertyOutput.
        ~Checkpointer();

        [[nodiscard]] bool ShouldWrite() const;

        void Write(configuration::SimConfig conf) const;

    private:
        void WriteConfig(configuration::SimConfig const&, std::filesystem::path const&) const;
        // Return the paths to the extraction and offset files
        [[nodiscard]] std::pair<path, path> WriteDists() const;
    };
}

#endif