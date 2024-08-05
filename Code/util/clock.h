// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UTIL_CLOCK_H
#define HEMELB_UTIL_CLOCK_H

namespace hemelb::util
{
    // Return a monotonic time in seconds from an arbitrary time point in the past
    double clock();
}
#endif
