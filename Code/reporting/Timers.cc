// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "reporting/Timers.hpp"

namespace hemelb::reporting
{
    // Explicit instantiations
    template class TimersBase<HemeLBClockPolicy>;
    template void TimersBase<HemeLBClockPolicy>::Reduce<net::IOCommunicator&>(net::IOCommunicator&);
}
