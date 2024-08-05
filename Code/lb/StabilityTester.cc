// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "lb/StabilityTester.h"

namespace hemelb::lb
{
    StabilityTester::StabilityTester(std::shared_ptr<const geometry::FieldData> iLatDat, net::Net* net,
                                     SimulationState* simState, reporting::Timers& timings,
                                     const hemelb::configuration::MonitoringConfig& testerConfig) :
            net::PhasedBroadcastRegular<>(net, simState, SPREADFACTOR), mLatDat(std::move(iLatDat)),
            mSimState(simState), timings(timings), testerConfig(testerConfig)
    {
        Reset();
    }

    bool StabilityTester::ShouldTerminateWhenConverged() const {
        return testerConfig.convergenceTerminate;
    }
    /**
     * Override the reset method in the base class, to reset the stability variables.
     */
    void StabilityTester::Reset()
    {
        mUpwardsStability = UndefinedStability;
        mDownwardsStability = UndefinedStability;

        mSimState->SetStability(UndefinedStability);

        std::fill(mChildrensStability.begin(), mChildrensStability.end(), UndefinedStability);
    }

    /**
     * Override the methods from the base class to propagate data from the root, and
     * to send data about this node and its childrens' stabilities up towards the root.
     */
    void StabilityTester::ProgressFromChildren(unsigned long splayNumber)
    {
        ReceiveFromChildren<int>(mChildrensStability.data(), 1);
    }

    void StabilityTester::ProgressFromParent(unsigned long splayNumber)
    {
        ReceiveFromParent<int>(&mDownwardsStability, 1);
    }

    void StabilityTester::ProgressToChildren(unsigned long splayNumber)
    {
        SendToChildren<int>(&mDownwardsStability, 1);
    }

    void StabilityTester::ProgressToParent(unsigned long splayNumber)
    {
        SendToParent<int>(&mUpwardsStability, 1);
    }

    /// Take the combined stability information (an int, with a value of hemelb::lb::Unstable
    /// if any child node is unstable) and start passing it back down the tree.
    void StabilityTester::TopNodeAction()
    {
        mDownwardsStability = mUpwardsStability;
    }

    void StabilityTester::PostReceiveFromChildren(unsigned long splayNumber)
    {
        timings.monitoring().Start();

        // No need to test children's stability if this node is already unstable.
        if (mUpwardsStability != Unstable)
        {
            if (std::any_of(
                    mChildrensStability.begin(), mChildrensStability.end(),
                    [](int _) {
                        return _ == Unstable;
                    }
            )) {
                mUpwardsStability = Unstable;
            }

            // If the simulation wasn't found to be unstable and we need to check for convergence, do it now.
            if ( (mUpwardsStability != Unstable) && testerConfig.doConvergenceCheck)
            {
                // mChildrensStability will contain UndefinedStability for non-existent children
                bool anyConverged = std::any_of(
                        mChildrensStability.begin(), mChildrensStability.end(),
                        [](int _) { return _ == StableAndConverged; }
                );
                bool anyStableNotConverged = std::any_of(
                        mChildrensStability.begin(), mChildrensStability.end(),
                        [](int _) { return _ == Stable; }
                );

                // With the current configuration the root node of the tree won't own any fluid sites. Its
                // state only depends on children nodes not on local state.
                if (anyConverged
                    && (mUpwardsStability == StableAndConverged || GetParent() == NOPARENT))
                {
                    mUpwardsStability = StableAndConverged;
                }

                if (anyStableNotConverged)
                {
                    mUpwardsStability = Stable;
                }
            }
        }

        timings.monitoring().Stop();
    }

    /// Apply the stability value sent by the root node to the simulation logic.
    void StabilityTester::Effect()
    {
        mSimState->SetStability((Stability) mDownwardsStability);
    }
}
