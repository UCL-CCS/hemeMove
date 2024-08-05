// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_IO_ENSURE_HEXFLOAT_H
#define HEMELB_IO_ENSURE_HEXFLOAT_H

#include <memory>

namespace hemelb::io {

    // If the current global locale doesn't support reading hexfloats,
    // add our one.
    //
    // Destructor restores the original.
    class GlobalHexFloatLocale {
        struct Impl;
        std::unique_ptr<Impl> impl;
    public:
        static bool CurrentCanParseHexFloats();
        GlobalHexFloatLocale();
        ~GlobalHexFloatLocale();
    };
}

#endif