// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_IO_FORMATS_RBC_H
#define HEMELB_IO_FORMATS_RBC_H

#include <cstddef>
#include <cstdint>
#include "io/formats/formats.h"

// Information on per-cell summary data
namespace hemelb::io::formats::rbc
{
    // Magic number to identify extraction data files.
    // ASCII for 'rbc' + EOF
    static constexpr std::uint32_t MagicNumber = 0x72626304;

    // The version number of the file format.
    static constexpr std::uint32_t VersionNumber = 1;

    // Header is 4 x uint32_t:
    // HemeLB magic, RBC magic, version, number of rows
    constexpr std::size_t header_size = 16U;

    // Each row is:
    //  - string representation of UUID (36 bytes)
    //  - cell barycentre (three doubles)
    constexpr std::size_t row_size = 36 + 3*8;
}
#endif
