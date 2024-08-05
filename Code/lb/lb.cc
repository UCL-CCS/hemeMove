// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "lb/lb.hpp"

namespace hemelb::lb {

    LBMBase::LBMBase(LbmParameters iParams, net::Net* net,
                     geometry::FieldData* latDat, SimulationState* simState,
                     reporting::Timers &atimings,
                     geometry::neighbouring::NeighbouringDataManager *neighbouringDataManager) :
            mNet(net), mLatDat(latDat), mState(simState),
            mParams(std::move(iParams)),
            timings(atimings), propertyCache(*simState, latDat->GetDomain()),
            neighbouringDataManager(neighbouringDataManager)
    {
    }

    void LBMBase::Initialise(BoundaryValues* iInletValues,
                                 BoundaryValues* iOutletValues)
    {
        mInletValues = iInletValues;
        mOutletValues = iOutletValues;

        InitCollisions();
    }

    void LBMBase::PrepareBoundaryObjects()
    {
        // First, iterate through all of the inlet and outlet objects, finding out the minimum density seen in the simulation.
        distribn_t minDensity = std::numeric_limits<distribn_t>::max();

        for (unsigned inlet = 0; inlet < mInletValues->GetLocalIoletCount(); ++inlet)
        {
            minDensity = std::min(minDensity, mInletValues->GetLocalIolet(inlet)->GetDensityMin());
        }

        for (unsigned outlet = 0; outlet < mOutletValues->GetLocalIoletCount(); ++outlet)
        {
            minDensity = std::min(minDensity, mOutletValues->GetLocalIolet(outlet)->GetDensityMin());
        }

        // Now go through them again, informing them of the minimum density.
        for (unsigned inlet = 0; inlet < mInletValues->GetLocalIoletCount(); ++inlet)
        {
            mInletValues->GetLocalIolet(inlet)->SetMinimumSimulationDensity(minDensity);
        }

        for (unsigned outlet = 0; outlet < mOutletValues->GetLocalIoletCount(); ++outlet)
        {
            mOutletValues->GetLocalIolet(outlet)->SetMinimumSimulationDensity(minDensity);
        }
    }

    void LBMBase::RequestComms()
    {
        timings.lb().Start();

        // Delegate to the lattice data object to post the asynchronous sends and receives
        // (via the Net object).
        // NOTE that this doesn't actually *perform* the sends and receives, it asks the Net
        // to include them in the ISends and IRecvs that happen later.
        mLatDat->SendAndReceive(mNet);

        timings.lb().Stop();
    }

    void LBMBase::EndIteration()
    {
    }
}