// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "tests/helpers/HasCommsTestFixture.h"

namespace hemelb::tests::helpers
{

    void HasCommsTestFixture::Init(const net::IOCommunicator& inst)
    {
        hemelbCommunicator = &inst;
    }

    const net::IOCommunicator& HasCommsTestFixture::Comms()
    {
        return *hemelbCommunicator;
    }

    const net::IOCommunicator* HasCommsTestFixture::hemelbCommunicator = nullptr;
}
