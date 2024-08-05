// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_HASSERT_H
#define HEMELB_HASSERT_H

#ifdef HEMELB_CODE
// A HemeLB assertion macro
#include "Exception.h"
#include "build_info.h"

#define HASSERT(expr) do {\
if constexpr (::hemelb::DEBUG) { \
    if (!(expr)) { \
        throw (::hemelb::Exception() << "Assertion failure '" #expr "' in '" __FILE__ ":" << __LINE__); \
    } \
}} while (0)

#else // Fall back to C assert macro
#include <cassert>
#define HASSERT(expr) assert(expr)

#endif

#endif
