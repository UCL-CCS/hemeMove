// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_TESTS_REPORTING_MOCKS_H
#define HEMELB_TESTS_REPORTING_MOCKS_H

#include <mpi.h>
#include "units.h"

namespace hemelb::net {
    class IOCommunicator;
}

namespace hemelb::tests
{
    class ClockMock {
    public:
        double operator()() {
            fakeTime += 10.0;
            return fakeTime;
        }

    private:
        double fakeTime = 0.0;
    };

    class MPICommsMock {
        unsigned calls = 1;
    public:
        inline int Size() const {
            return 5;
        }

        template <std::size_t N = std::dynamic_extent>
        void Reduce(std::span<double, N> recvbuf, std::span<double const, N> sendbuf, MPI_Op, int root) {
            auto const count = sendbuf.size();
            //REQUIRE(reporting::Timers::numberOfTimers == count);

            for (int i = 0; i < count; i++) {
                REQUIRE(10.0 * i == sendbuf[i]);
                recvbuf[i] = 5.0 * i * calls;
            }
            calls++;
        }
    };
}

#endif // HEMELB_TESTS_REPORTING_MOCKS_H
