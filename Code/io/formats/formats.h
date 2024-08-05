// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_IO_FORMATS_FORMATS_H
#define HEMELB_IO_FORMATS_FORMATS_H

#include <cstdint>

namespace hemelb::io::formats
{
    // This is the generic HemeLB binary file identifier. It should be the
    // first 4 bytes of every (binary) file used for IO. It should be then
    // followed by another 4 bytes identifying the particular type/version,
    // that number being terminated by the EOF character (0x04)
    constexpr std::uint32_t HemeLbMagicNumber = 0x686c6221;
}
#endif
