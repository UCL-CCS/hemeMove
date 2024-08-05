// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_IO_TIMEPATTERN_H
#define HEMELB_IO_TIMEPATTERN_H

#include <string>
#include <string_view>

#include "units.h"

namespace hemelb::io {
    // We have, for output paths the concept of a string or path that needs
    // to be templated on the time step.
    // In the XML config file, this should be represented as "%d".

    class TimePattern {
        std::string pattern;
        std::size_t iPercent = 0;

    public:
        // Check that the string is OK
        static bool Check(std::string_view p);

        TimePattern() = default;
        TimePattern(std::string_view p);

        // Remove the time pattern part
        std::string Strip() const;
        std::string Format(LatticeTimeStep t, LatticeTimeStep max) const;
    };
}
#endif
