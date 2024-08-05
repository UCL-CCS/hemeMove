// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_ENTROPYTESTER_H
#define HEMELB_LB_ENTROPYTESTER_H

#include "net/PhasedBroadcastRegular.h"
#include "geometry/Domain.h"
#include "lb/HFunction.h"
#include "log/Logger.h"

namespace hemelb::lb
{
    template<lattice_type LatticeType>
    class EntropyTester : public net::PhasedBroadcastRegular<false, 1, 1, false, true>
    {
    private:
        enum HTHEOREM
        {
            OBEYED,
            DISOBEYED
        };

        /**
         * Slightly arbitrary spread factor for the tree.
         */
        static constexpr unsigned SPREADFACTOR = 10;

        geometry::Domain const* mLatDat;

        /**
         * Stability value of this node and its children to propagate upwards.
         */
        int mUpwardsValue = OBEYED;

        static constexpr auto ObeyedArray() {
            std::array<int, SPREADFACTOR> x;
            std::fill(x.begin(), x.end(), OBEYED);
            return x;
        }

        /**
         * Array for storing the passed-up stability values from child nodes.
         */
        std::array<int, SPREADFACTOR> mChildrensValues = ObeyedArray();

        std::array<bool, COLLISION_TYPES> mCollisionTypesTested;
        std::unique_ptr<double[]> mHPreCollision;

    public:
        EntropyTester(int* collisionTypes, unsigned int typesTested,
                      const geometry::Domain * iLatDat, net::Net* net,
                      SimulationState* simState) :
                net::PhasedBroadcastRegular<false, 1, 1, false, true>(net, simState, SPREADFACTOR),
                mLatDat(iLatDat)
        {
            std::fill(mCollisionTypesTested.begin(), mCollisionTypesTested.end(), false);
            for (unsigned int i = 0; i < typesTested; i++)
            {
                mCollisionTypesTested[collisionTypes[i]] = true;
            }

            mHPreCollision.reset(new double[mLatDat->GetLocalFluidSiteCount()]);

            Reset();
        }

        ~EntropyTester() override = default;

        void PreReceive() override
        {
            double dHMax = 0.0;

            // The order of arguments in max is important
            // If distributions go negative HFunc.eval() will return a NaN. Due to the
            // nature of NaN and the structure of max, dH will be assigned as NaN if this is the case
            // This is what we want, because the EntropyTester will fail otherwise and abort when
            // it is simply sufficient to wait until StabilityTester restarts.

            for (unsigned collision_type = 0; collision_type < COLLISION_TYPES; collision_type++) {
                if (mCollisionTypesTested[collision_type]) {
                    auto [begin, end] = mLatDat->GetMidDomainSiteRange(collision_type);
                    for (site_t i = begin; i < end; ++i) {
                        auto const site = mLatDat->GetSite(i);

                        HFunction<LatticeType> HFunc(site.GetFOld<LatticeType>(), nullptr);
                        dHMax = std::max(dHMax, HFunc.eval() - mHPreCollision[i]);
                    }
                }
            }

            for (unsigned collision_type = 0; collision_type < COLLISION_TYPES; collision_type++) {
                auto [begin, end] = mLatDat->GetDomainEdgeSiteRange(collision_type);
                if (mCollisionTypesTested[collision_type]) {
                    for (site_t i = begin; i < end; ++i) {
                        auto const site = mLatDat->GetSite(i);

                        HFunction<LatticeType> HFunc(site.GetFOld<LatticeType>(), nullptr);
                        dHMax = std::max(dHMax, HFunc.eval() - mHPreCollision[i]);
                    }
                }
            }

            // Ideally dH should never be greater than zero. However, because H is positive
            // definite accuracy as well as rounding and truncation errors can make dH greater
            // than zero in certain cases. The tolerance is limited by the accuracy of the
            // simulation (including the accuracy to which alpha is calculated) The tolerance
            // has to be at least as big as the accuracy to which alpha is calculated.
            if (dHMax > 1.0E-6)
                mUpwardsValue = DISOBEYED;
        }

        // Override the reset method in the base class, to reset the stability variables.
        void Reset()
        {
            // Re-initialise all values to indicate that the H-theorem is obeyed.
            mUpwardsValue = OBEYED;
            mChildrensValues = ObeyedArray();
        }

    protected:
        // Override the methods from the base class to propagate data from the root, and
        // to send data about stability of this node and its children up towards the root.
        void ProgressFromChildren(unsigned long splayNumber) override
        {
            ReceiveFromChildren<int>(mChildrensValues.data(), 1);
        }

        void ProgressToParent(unsigned long splayNumber) override
        {
            for (unsigned collision_type = 0; collision_type < COLLISION_TYPES; collision_type++) {
                if (mCollisionTypesTested[collision_type]) {
                    auto [begin, end] = mLatDat->GetMidDomainSiteRange(collision_type);
                    for (site_t i = begin; i < end; ++i) {
                        auto const site = mLatDat->GetSite(i);
                        HFunction<LatticeType> HFunc(site.GetFOld<LatticeType>(), nullptr);
                        mHPreCollision[i] = HFunc.eval();
                    }
                }
            }

            for (unsigned collision_type = 0; collision_type < COLLISION_TYPES; collision_type++) {
                if (mCollisionTypesTested[collision_type]) {
                    auto [begin, end] = mLatDat->GetDomainEdgeSiteRange(collision_type);
                    for (site_t i = begin; i < end; ++i) {
                        auto const site = mLatDat->GetSite(i);
                        HFunction<LatticeType> HFunc(site.GetFOld<LatticeType>(), nullptr);
                        mHPreCollision[i] = HFunc.eval();
                    }
                }
            }
            SendToParent<int>(&mUpwardsValue, 1);
        }

        /**
         * Take the combined stability information (an int, with a value of hemelb::lb::Unstable
         * if any child node is unstable) and start passing it back down the tree.
         */
        void TopNodeAction() override
        {
            if (mUpwardsValue == DISOBEYED)
            {
                log::Logger::Log<log::Error, log::Singleton>("H Theorem violated.");
            }
        }

        // Override the method from the base class to use the data from child nodes.
        void PostReceiveFromChildren(unsigned long splayNumber) override {
            // No need to test children's entropy direction if this node already disobeys H-theorem.
            if (mUpwardsValue == OBEYED) {
                if (std::any_of(
                        mChildrensValues.begin(), mChildrensValues.end(),
                        [](int _) { return _ == DISOBEYED; }
                )) {
                    mUpwardsValue = DISOBEYED;
                }
            }
        }
    };
}

#endif /* HEMELB_LB_ENTROPYTESTER_H */
