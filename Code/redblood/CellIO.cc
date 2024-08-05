// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "redblood/CellIO.h"

#include <concepts>
#include <cstring>
#include <boost/uuid/uuid_io.hpp>

#include "configuration/PathManager.h"
#include "io/formats/rbc.h"
#include "io/writers/XdrWriter.h"
#include "io/readers/XdrFileReader.h"
#include "io/readers/XdrMemReader.h"
#include "lb/SimulationState.h"
#include "log/Logger.h"
#include "net/MpiFile.h"
#include "redblood/CellBase.h"
#include "redblood/FaderCell.h"
#include "redblood/MeshIO.h"
#include "util/Iterator.h"
#include "util/span.h"

namespace hemelb::redblood {

    namespace {
        static_assert(std::convertible_to<CellVtkOutput, CellChangeListener>);

        // Create output directory for current writing step. Requires syncing to
        // ensure no process goes ahead before directory is created.
        [[nodiscard]] std::filesystem::path EnsureOutputDir(CellOutputBase const &outputter) {
            return outputter.fileManager->GetRBCOutputPathWithSubdir(
                    std::to_string(outputter.simState->GetTimeStep()),
                    true,
                    true
            );
        }

        constexpr std::size_t UUID_STRING_LEN = 36;
    }

    void CellVtkOutput::operator()(CellContainer const &cells) {
        auto timestep = simState->GetTimeStep();
        if ((timestep % period) != 0)
            return;

        log::Logger::Log<log::Info, log::OnePerCore>("printstep %d, num cells %d", timestep, cells.size());

        auto rbcOutputDir = EnsureOutputDir(*this);
        ioComms.Barrier();

        for (auto &cell: cells) {
            // to_chars guarantees to write exactly 36 chars, also need .vtp\0
            char name[UUID_STRING_LEN + 5];
            boost::uuids::to_chars(cell->GetTag(), name);
            std::strncpy(name + UUID_STRING_LEN, ".vtp", 5);
            auto filename = rbcOutputDir / name;
            std::shared_ptr<redblood::CellBase> cell_base = [&cell]() {
                if (auto fader = std::dynamic_pointer_cast<redblood::FaderCell>(cell)) {
                    return fader->GetWrapeeCell();
                } else {
                    return cell;
                }
            }();

            auto cell_cast = std::dynamic_pointer_cast<redblood::Cell>(cell_base);
            assert(cell_cast);
            auto meshio = redblood::VTKMeshIO{};
            meshio.writeFile(filename.native(), *cell_cast, unitConverter.get());
        }
    }

    CellBarycentreOutput::CellBarycentreOutput(hemelb::LatticeTimeStep period,
                                               std::shared_ptr<const util::UnitConverter> unitConverter,
                                               std::shared_ptr<const lb::SimulationState> simState,
                                               std::shared_ptr<const configuration::PathManager> fileManager,
                                               net::IOCommunicator comms) :
            CellOutputBase{
                    period, std::move(unitConverter), std::move(simState), std::move(fileManager), comms
    } {
    }

    void CellBarycentreOutput::operator()(CellContainer const& cells) {
        namespace fmt = io::formats;
        auto timestep = simState->GetTimeStep();
        if ((timestep % period) != 0)
            return;

        auto rbcOutputDir = EnsureOutputDir(*this);
        // Barrier before use!

        unsigned n_local_cells = cells.size();
        unsigned n_cells_before_me = 0U;
        auto req = ioComms.Iexscan(n_local_cells, n_cells_before_me, MPI_SUM);
        const std::size_t local_data_start = ioComms.OnIORank() ? fmt::rbc::header_size : 0;
        const std::size_t local_write_size =  local_data_start + n_local_cells * fmt::rbc::row_size;
        std::vector<std::byte> buffer(local_write_size);

        // Build this rank's part of the main file,
        // will fill header later once know the total number of cells
        auto xdrWriter = io::MakeXdrWriter(buffer.begin() + local_data_start, buffer.end());

        for (auto [i, cell]: util::enumerate(cells)) {
            std::byte tag[UUID_STRING_LEN];
            boost::uuids::to_chars(cell->GetTag(), reinterpret_cast<char*>(tag));
            xdrWriter << std::span<std::byte>(tag, UUID_STRING_LEN) << cell->GetBarycentre();
        }

        auto bary_filename = rbcOutputDir / "barycentres.rbc";

        req.Wait();

        // On last rank, know the total number of cells. Pass to rank 0 (unless we are rank 0 too!)
        int const last_rank = ioComms.Size() - 1;
        auto total_cells = n_cells_before_me + n_local_cells;
        if (last_rank != 0) {
            if (ioComms.Rank() == last_rank) {
                ioComms.Send(total_cells, 0);
            }
            if (ioComms.Rank() == 0) {
                ioComms.Receive(total_cells, last_rank);
            }
        }

        if (ioComms.OnIORank()) {
            // Create header
            auto headerWriter = io::MakeXdrWriter(buffer.begin(), buffer.begin() + fmt::rbc::header_size);
            headerWriter << fmt::HemeLbMagicNumber << fmt::rbc::MagicNumber << fmt::rbc::VersionNumber
                         << std::uint32_t(total_cells);
        }

        // Now know our write position and that directory is created.
        const std::size_t local_write_start = ioComms.OnIORank() ? 0 : (fmt::rbc::header_size + n_cells_before_me * fmt::rbc::row_size);

        auto bary_file = net::MpiFile::Open(ioComms, bary_filename, MPI_MODE_WRONLY | MPI_MODE_CREATE | MPI_MODE_EXCL);
        bary_file.WriteAt(local_write_start, to_const_span(buffer));
        bary_file.Close();
    }

    std::uint32_t CellBarycentreInput::ReadHeader() const {
        namespace fmt = io::formats;
        io::XdrFileReader headerReader(filename);
        std::uint32_t hlbMagicNumber, rbcMagicNumber, version, nCells;
        headerReader.read(hlbMagicNumber);
        headerReader.read(rbcMagicNumber);
        headerReader.read(version);
        headerReader.read(nCells);

        if (hlbMagicNumber != fmt::HemeLbMagicNumber)
            throw Exception() << "This file does not start with the HemeLB magic number."
                              << " Expected: " << unsigned(fmt::HemeLbMagicNumber)
                              << " Actual: " << hlbMagicNumber;

        if (rbcMagicNumber != fmt::rbc::MagicNumber)
            throw Exception() << "This file does not have the offset magic number."
                              << " Expected: " << unsigned(fmt::rbc::MagicNumber)
                              << " Actual: " << rbcMagicNumber;

        if (version != fmt::rbc::VersionNumber)
            throw Exception() << "Version number incorrect."
                              << " Supported: " << unsigned(fmt::rbc::VersionNumber)
                              << " Input: " << version;
        return nCells;
    }

    auto CellBarycentreInput::ReadRows(net::MpiCommunicator const& comms, std::size_t start, std::size_t end) const -> std::vector<row> {
        namespace fmt = io::formats;
        auto const n = end - start;
        auto inputFile = net::MpiFile::Open(comms, filename, MPI_MODE_RDONLY);
        auto read_start = fmt::rbc::header_size + start * fmt::rbc::row_size;
        auto read_size = n * fmt::rbc::row_size;

        std::vector<std::byte> read_buf(read_size);
        inputFile.ReadAt(read_start, to_span(read_buf));
        auto reader = io::XdrMemReader(read_buf);
        std::vector<row> ans(n);
        boost::uuids::string_generator gen;
        for (auto i = 0; i < n; ++i) {
            std::byte tag[UUID_STRING_LEN];
            auto s = std::span{tag};
            reader.read(s);
            ans[i].first = gen(reinterpret_cast<char*>(tag), reinterpret_cast<char*>(tag)+UUID_STRING_LEN);
            reader.read(ans[i].second);
        }
        return ans;
    }
}