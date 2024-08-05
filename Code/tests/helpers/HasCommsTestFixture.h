// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_TESTS_HELPERS_HASCOMMSTESTFIXTURE_H
#define HEMELB_TESTS_HELPERS_HASCOMMSTESTFIXTURE_H

#include "net/IOCommunicator.h"

namespace hemelb::tests::helpers
{
    class HasCommsTestFixture
    {
    private:
        static net::IOCommunicator const* hemelbCommunicator;
    public:
        static void Init(const net::IOCommunicator& inst);
    protected:
        static const net::IOCommunicator& Comms();
    };
}
#endif // HEMELB_TESTS_HELPERS_HASCOMMSTESTFIXTURE_H
