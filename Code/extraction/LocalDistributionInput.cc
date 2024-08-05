// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#include "extraction/LocalDistributionInput.h"

#include "extraction/OutputField.h"
#include "geometry/FieldData.h"
#include "io/formats/formats.h"
#include "io/formats/extraction.h"
#include "io/formats/offset.h"
#include "io/readers/XdrFileReader.h"
#include "io/readers/XdrMemReader.h"
#include "log/Logger.h"
#include "util/span.h"

namespace hemelb::extraction {
    namespace fmt = hemelb::io::formats;

    LocalDistributionInput::LocalDistributionInput(std::filesystem::path dataFilePath,
                                                   std::optional<std::filesystem::path> maybeOffsetPath,
                                                   const net::IOCommunicator& ioComm) :
            comms{ioComm}, filePath{std::move(dataFilePath)}
    {
        if (maybeOffsetPath) {
            offsetPath = std::move(*maybeOffsetPath);
        } else {
            offsetPath = fmt::offset::ExtractionToOffset(filePath);
        }
    }

    namespace {
        // The required xtr field header len
        uint64_t constexpr expectedFieldHeaderLength = 40U;
        uint64_t constexpr totalXtrHeaderLength = fmt::extraction::MainHeaderLength + expectedFieldHeaderLength;
    }

    void LocalDistributionInput::LoadDistribution(geometry::FieldData* latDat, std::optional<LatticeTimeStep>& targetTime)
    {
        auto&& dom = latDat->GetDomain();
        const auto NUMVECTORS = dom.GetLatticeInfo().GetNumVectors();

        // We could supply hints regarding how the file should be read
        // but we are not doing so yet.

        // Open the file as read-only.
        // TODO: raise an exception if the file does not exist.
        auto inputFile = net::MpiFile::Open(comms, filePath, MPI_MODE_RDONLY);
        // Set the view to the file.
        inputFile.SetView(0, MPI_CHAR, MPI_CHAR, "native");
        ReadExtractionHeaders(inputFile, NUMVECTORS);

        // Now read offset file.
        ReadOffsets(offsetPath);

        // Figure out how many checkpoints are in the XTR file and
        // therefore the position to start at.
        auto nTimes = [&](){
            uint64_t fileSize = inputFile.GetSize();
            auto dataSize = fileSize - totalXtrHeaderLength;
            if (dataSize % allCoresWriteLength)
                throw Exception() << "Checkpoint file length not consistent with integer number of checkpoints";
            return dataSize / allCoresWriteLength;
        }();

        auto ReadTimeByIndex = [&](uint64_t iTS) {
            uint64_t ans;
            std::vector<std::byte> tsbuf(8);
            inputFile.ReadAt(localStart + iTS*allCoresWriteLength, to_span(tsbuf));
            io::XdrMemReader dataReader(tsbuf);
            dataReader.read(ans);
            return ans;
        };

        uint64_t iTS;
        if (comms.OnIORank()) {
            if (targetTime) {
                // We have a target time - look for it in the file
                iTS = 0;
                uint64_t len = nTimes;
                while(len != 0) {
                    auto l2 = len/2;
                    auto m = iTS + l2;
                    timestep = ReadTimeByIndex(m);
                    if (timestep < *targetTime) {
                        iTS = m + 1;
                        len -= l2 + 1;
                    } else {
                        len = l2;
                    }
                }

                if (timestep != *targetTime)
                    throw Exception() << "Target timestep " << *targetTime << " not found in checkpoint file.";
            } else {
                // initial time unspecified, use the last one
                iTS = nTimes - 1;
                timestep = ReadTimeByIndex(iTS);
            }
        }
        comms.Broadcast(timestep, comms.GetIORank());
        comms.Broadcast(iTS, comms.GetIORank());
        if (!targetTime)
            targetTime = timestep;

        log::Logger::Log<log::Info, log::Singleton>("Reading checkpoint from timestep %d with index %d", timestep, iTS);
        // Read the local part of the checkpoint
        const auto readLength = localStop - localStart;
        const auto timeStart = iTS * allCoresWriteLength;

        std::vector<std::byte> dataBuffer(readLength);
        inputFile.ReadAt(timeStart + localStart, to_span(dataBuffer));
        io::XdrMemReader dataReader(dataBuffer);

        // Read the timestep
        if (comms.OnIORank()) {
            dataReader.read(timestep);
        }

        site_t iSite = 0;
        for (; dataReader.GetPosition() < readLength; iSite++) {
            // Read the grid coord and check it's consistent with latDat.
            {
                // Stored as 32 b unsigned
                util::Vector3D<uint32_t> tmp;
                dataReader.read(tmp.x());
                dataReader.read(tmp.y());
                dataReader.read(tmp.z());

                // Convert to canonical type
                util::Vector3D<site_t> grid{tmp};
                // Look up the site ID and rank, as decomposed by this run
                // of HemeLB, for the grid coordinate read from the
                // checkpoint file.
                proc_t rank; site_t index;
                if (!dom.GetContiguousSiteId(grid, rank, index)) {
                    // function returns a 'valid' flag
                    throw Exception() << "Cannot get valid site from extracted site coordinate";
                }
                if (rank != comms.Rank())
                    throw Exception() << "Site read on rank " << comms.Rank()
                                      << " but should be read on " << rank;
                if (index != iSite)
                    throw Exception() << "Site read at index " << iSite
                                      << " but should be read at " << index;
            }

            distribn_t* f_old_p = latDat->GetFOld(iSite * NUMVECTORS);
            distribn_t* f_new_p = latDat->GetFNew(iSite * NUMVECTORS);
            // distField is read on IO rank and checked to be equal to
            // NUMVECTORS so we use that instead of broadcasting and
            // storing.
            for (auto i = 0U; i < NUMVECTORS; i++) {
                distribn_t field_val;
                dataReader.read(field_val);
                f_new_p[i] = f_old_p[i] = field_val;
            }
        }

        if (iSite != dom.GetLocalFluidSiteCount())
            throw Exception() << "Read " << iSite
                              << " sites but expected " << dom.GetLocalFluidSiteCount();
    }

    void LocalDistributionInput::ReadExtractionHeaders(net::MpiFile& inputFile, const unsigned NUMVECTORS) {
        // The headers technically aren't needed (because of the offset
        // file), but we check that they are as expected.
        if (comms.OnIORank()) {
            auto preambleBuf = std::vector<std::byte>(fmt::extraction::MainHeaderLength);
            inputFile.Read(to_span(preambleBuf));
            auto preambleReader = io::XdrMemReader(preambleBuf);

            // Read the magic numbers.
            uint32_t hlbMagicNumber, extMagicNumber, version;
            preambleReader.read(hlbMagicNumber);
            preambleReader.read(extMagicNumber);
            preambleReader.read(version);

            // Check the value of the HemeLB magic number.
            if (hlbMagicNumber != fmt::HemeLbMagicNumber)
            {
                throw Exception() << "This file does not start with the HemeLB magic number."
                                  << " Expected: " << unsigned(fmt::HemeLbMagicNumber)
                                  << " Actual: " << hlbMagicNumber;
            }

            // Check the value of the extraction file magic number.
            if (extMagicNumber != fmt::extraction::MagicNumber)
            {
                throw Exception() << "This file does not have the extraction magic number."
                                  << " Expected: " << unsigned(fmt::extraction::MagicNumber)
                                  << " Actual: " << extMagicNumber;
            }

            // Check the version number.
            if (version != fmt::extraction::VersionNumber)
            {
                throw Exception() << "Version number incorrect."
                                  << " Supported: " << unsigned(fmt::extraction::VersionNumber)
                                  << " Input: " << version;
            }

            // Perhaps should check these...
            {
                // Obtain the size of voxel in metres.
                double dx, dt, dm;
                preambleReader.read(dx);
                preambleReader.read(dt);
                preambleReader.read(dm);

                // Obtain the origin.
                double origin[3];
                preambleReader.read(origin[0]);
                preambleReader.read(origin[1]);
                preambleReader.read(origin[2]);
                double p0;
                preambleReader.read(p0);
            }
            // Obtain the total number of sites, fields & header len
            uint64_t numberOfSites;
            uint32_t numberOfFields, lengthOfFieldHeader;
            preambleReader.read(numberOfSites);
            preambleReader.read(numberOfFields);
            preambleReader.read(lengthOfFieldHeader);

            if (numberOfFields != 1 )
                throw Exception() << "Checkpoint file must contain exactly one field, the distributions, but has "
                                  << numberOfFields;
            if (lengthOfFieldHeader != expectedFieldHeaderLength)
                throw Exception() << "Checkpoint file's field header must be "
                                  << expectedFieldHeaderLength << " B long, but is "
                                  << lengthOfFieldHeader << " B";

            auto fieldHeaderBuf = std::vector<std::byte>(lengthOfFieldHeader);
            inputFile.Read(to_span(fieldHeaderBuf));
            auto fieldHeaderReader = io::XdrMemReader(fieldHeaderBuf);

            fieldHeaderReader.read(distField.name);
            if (distField.name != "distributions")
                throw Exception() << "Checkpoint file must contain field named 'distributions', but has '"
                                  << distField.name << "'";

            fieldHeaderReader.read(distField.numberOfElements);
            if (distField.numberOfElements != NUMVECTORS)
                throw Exception() << "Checkpoint field distributions contains " << distField.numberOfElements
                                  << " distributions but this build of HemeLB requires " << NUMVECTORS;

            fieldHeaderReader.read(distField.typecode);
            if (distField.typecode != static_cast<std::uint32_t>(io::formats::extraction::TypeCode::DOUBLE))
                throw Exception() << "Checkpoint contains wrong data type";

            fieldHeaderReader.read(distField.numberOfOffsets);
            if (distField.numberOfOffsets != 0)
                throw Exception() << "Checkpoint should not have offsets";

            double scale;
            fieldHeaderReader.read(scale);
            if (scale != 0.0)
                throw Exception() << "Checkpoint has scaling applied!";

        }
    }

    void LocalDistributionInput::ReadOffsets(const std::string& offsetFileName) {
        std::vector<uint64_t> offsets;

        // Only actually read on IO rank
        if (comms.OnIORank()) {
            io::XdrFileReader offsetReader(offsetFileName);
            uint32_t hlbMagicNumber, offMagicNumber, version;
            int32_t nRanks;
            offsetReader.read(hlbMagicNumber);
            offsetReader.read(offMagicNumber);
            offsetReader.read(version);
            offsetReader.read(nRanks);

            if (hlbMagicNumber != fmt::HemeLbMagicNumber)
                throw Exception() << "This file does not start with the HemeLB magic number."
                                  << " Expected: " << unsigned(fmt::HemeLbMagicNumber)
                                  << " Actual: " << hlbMagicNumber;

            if (offMagicNumber != fmt::offset::MagicNumber)
                throw Exception() << "This file does not have the offset magic number."
                                  << " Expected: " << unsigned(fmt::offset::MagicNumber)
                                  << " Actual: " << offMagicNumber;

            if (version != fmt::offset::VersionNumber)
                throw Exception() << "Version number incorrect."
                                  << " Supported: " << unsigned(fmt::offset::VersionNumber)
                                  << " Input: " << version;

            if (nRanks != comms.Size())
                throw Exception() << "Offset file has wrong number of MPI ranks."
                                  << " Running with: " << comms.Size()
                                  << " Input: " << nRanks;

            // Now read the encoded nProcs+1 values
            // We are going to duplicate these into a flattened array of shape (nRanks, 2)
            // [start0, end0, start1, end1, ...]
            // where end_i == start_i+1 (except for the start finish obvs)
            offsets.resize(2*nRanks);
            offsetReader.read(offsets[0]);
            for (int i = 1; i < nRanks; ++i) {
                offsetReader.read(offsets[2*i]);
                offsets[2*i - 1] = offsets[2*i];
            }
            offsetReader.read(offsets[2*nRanks-1]);
            // Compute the total length of a record
            allCoresWriteLength = offsets[2*nRanks-1] - offsets[0];
        }

        // Now bcast/scatter from IO rank to all
        comms.Broadcast(allCoresWriteLength, comms.GetIORank());
        auto start_finish = comms.Scatter(offsets, 2, comms.GetIORank());
        localStart = start_finish[0];
        localStop = start_finish[1];
    }
}
