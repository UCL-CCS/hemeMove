// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_NET_PHASED_CONCERN_H
#define HEMELB_NET_PHASED_CONCERN_H
namespace hemelb::net::phased
{
    class Concern
    {
    public:
        Concern() = default;
        virtual ~Concern() = default;
        virtual bool CallAction(int action) = 0;
    };
}
#endif
