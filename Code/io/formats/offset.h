// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_IO_FORMATS_OFFSET_H
#define HEMELB_IO_FORMATS_OFFSET_H

namespace hemelb::io::formats
{
    struct offset {
        // Magic number to identify offset files.
        // ASCII for 'off' + EOF
        static constexpr std::uint32_t MagicNumber = 0x6f666604;

        // The version number of the file format.
        static constexpr std::uint32_t VersionNumber = 1;

        // Header contains:
        // - HemeLb magic - uint32
        // - Offset magic - uint32
        // - Offset version - uint32
        // - number of ranks - uint32
        static constexpr std::size_t HeaderLength = 16;

        // Body contains, an array of size (number of ranks + 1),
        // where elem[i] contains the offset for rank i and elem[i+1]
        // contains the past the end element for that rank.
        static constexpr std::size_t RecordLength = sizeof(uint64_t);

        // Helper function
        static inline std::string ExtractionToOffset(const std::string &xtrPath) {
            auto iDot = xtrPath.rfind('.');
            if (iDot == std::string::npos)
                throw Exception() << "Cannot split extension from extraction filename";
            return xtrPath.substr(0, iDot) + ".off";
        }
    };
}
#endif
