// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_REPORTING_POLICIES_H
#define HEMELB_REPORTING_POLICIES_H

/**
 * @file
 * This file contains Policy classes defining how the report generator should interact with the file system, MPI, and the clock.
 * These are used as template arguments to the principal classes.
 * In file unittests/reporting/Mocks.h, mock versions of these policies are defined to facilitate testing.
 */

#include <fstream>
#include "net/mpi.h"
#include "net/IOCommunicator.h"
#include "util/numerical.h"

namespace hemelb::reporting
{
    /**
     * A way to get the time.
     * Mocked by hemelb::unittests::reporting::ClockMock
     */
    struct HemeLBClockPolicy
    {
        /**
         * Get the time
         * @return current time in seconds.
         */
        inline double operator()() const
        {
            return MPI_Wtime();
        }
    };
    static_assert(std::same_as<std::invoke_result_t<HemeLBClockPolicy>, double>);
}
#endif
