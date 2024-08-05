// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "io/TimePattern.h"

#include <string>

#include "Exception.h"

namespace hemelb::io {
    bool TimePattern::Check(std::string_view p) {
        // Path must contain exactly one printf conversion specifier
        // for an integer.
        auto n_pc = std::count(p.begin(), p.end(), '%');
        auto i_pcd = p.find("%d", 0, 2);
        return (n_pc == 1) && (i_pcd != std::string::npos);
    }

    TimePattern::TimePattern(std::string_view p) : pattern() {
        if (!Check(p))
            throw (Exception() << "Invalid pattern");
        iPercent = p.find("%d", 0, 2);
        // The part before %d
        auto beginning = p.substr(0, iPercent);
        // The part after
        auto end = p.substr(iPercent + 2);
        // Build the pattern
        pattern.reserve(p.size() + 2);
        pattern += beginning;
        pattern += "%*ld";
        pattern += end;
    }

    std::string TimePattern::Strip() const {
        std::string ans;
        ans.reserve(pattern.size() - 4);
        ans += pattern.substr(0, iPercent);
        ans += pattern.substr(iPercent + 4);
        return ans;
    }

    std::string TimePattern::Format(LatticeTimeStep t, LatticeTimeStep max) const {
        int prec = 3;
        unsigned long next = 1000;
        while (max > next) {
            prec += 1;
            next *= 10;
        }

        int sz = std::snprintf(nullptr, 0,
                               pattern.data(), prec, t);
        if (sz < 0)
            throw (Exception() << "Formatting error");

        // +1 for the null terminator
        std::string ans(sz, '\0');
        std::snprintf(ans.data(), ans.size() + 1,
                      pattern.data(), prec, t);
        return ans;
    }
}