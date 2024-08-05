// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_REDBLOOD_CELLIO_H
#define HEMELB_REDBLOOD_CELLIO_H

#include <array>
#include <filesystem>
#include <memory>
#include <boost/uuid/uuid.hpp>

#include "units.h"
#include "net/IOCommunicator.h"
#include "redblood/types_fwd.h"

namespace hemelb::configuration { class PathManager; }
namespace hemelb::lb { class SimulationState; }
namespace hemelb::util { class UnitConverter; }

namespace hemelb::redblood {

    // Base for outputting RBC data.
    //
    // Derived types must implement operator()(CellContainer const&)
    struct CellOutputBase {
        // How ofter to write
        LatticeTimeStep period;
        // nullptr => use lattice units
        std::shared_ptr<util::UnitConverter const> unitConverter;
        std::shared_ptr<lb::SimulationState const> simState;
        std::shared_ptr<configuration::PathManager const> fileManager;
        net::IOCommunicator ioComms;
    };

    // Write cells as vtkPolyData (.vtp)
    struct CellVtkOutput : CellOutputBase {
        void operator()(CellContainer const& cells);
    };

    // Write cell barycentres in HemeLB binary (.rbc)
    class CellBarycentreOutput : CellOutputBase {
    public:
        CellBarycentreOutput(
                LatticeTimeStep period,
                std::shared_ptr<util::UnitConverter const> unitConverter,
                std::shared_ptr<lb::SimulationState const> simState,
                std::shared_ptr<configuration::PathManager const> fileManager,
                net::IOCommunicator comms
        );

        void operator()(CellContainer const& cells);
    };

    struct CellBarycentreInput {
        std::filesystem::path filename;
        [[nodiscard]] std::uint32_t ReadHeader() const;
        using uuid = boost::uuids::uuid;
        using row = std::pair<uuid, LatticePosition>;
        [[nodiscard]] std::vector<row> ReadRows(net::MpiCommunicator const& comms, std::size_t start, std::size_t end) const;
    };
}

#endif
