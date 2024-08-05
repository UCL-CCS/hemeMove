// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_LB_H
#define HEMELB_LB_LB_H

#include "net/net.h"
#include "net/IteratedAction.h"
#include "net/IOCommunicator.h"
#include "lb/SimulationState.h"
#include "lb/InitialCondition.h"
#include "lb/iolets/BoundaryValues.h"
#include "lb/MacroscopicPropertyCache.h"
#include "util/UnitConverter.h"
#include "reporting/Timers.h"
#include "Traits.h"
#include <typeinfo>

/**
 * Namespace 'lb' contains classes for the scientific core of the Lattice Boltzmann simulation
 */
namespace hemelb::lb
{

    /// Holds LBM data and behaviour that does not depend on the velocity set/collision etc.
    class LBMBase : public net::IteratedAction {
    protected:
        net::Net* mNet;
        geometry::FieldData* mLatDat;
        SimulationState* mState;
        BoundaryValues *mInletValues = nullptr;
        BoundaryValues* mOutletValues = nullptr;

        LbmParameters mParams;

        hemelb::reporting::Timers &timings;

        MacroscopicPropertyCache propertyCache;

        geometry::neighbouring::NeighbouringDataManager *neighbouringDataManager;

    public:
        /// OK, with apologies we have a 2-stage construction process.
        /// The result of this constructor isn't yet ready to simulate, it must
        /// have Initialise(...) called also.
        ///
        /// Constructor separated due to need to access the partially initialized
        /// LBM in order to initialize the arguments to the second construction phase.
        LBMBase(LbmParameters params, net::Net* net,
            geometry::FieldData* latDat, SimulationState* simState, reporting::Timers &atimings,
            geometry::neighbouring::NeighbouringDataManager *neighbouringDataManager);

        /// Second constructor.
        void Initialise(BoundaryValues* iInletValues,
                        BoundaryValues* iOutletValues);

        void RequestComms() override; ///< part of IteratedAction interface.
        void EndIteration() override; ///< part of IteratedAction interface.

        inline LbmParameters* GetLbmParams() {
            return &mParams;
        }
        inline MacroscopicPropertyCache& GetPropertyCache() {
            return propertyCache;
        }

        /// Interface to apply initial conditions
        virtual void SetInitialConditions(lb::InitialCondition const& ic_conf, const net::IOCommunicator& ioComms) = 0;

    protected:
        virtual void InitCollisions() = 0;

        /**
         * Ensure that the BoundaryValues objects have all necessary fields populated.
         */
        void PrepareBoundaryObjects();
    };

    /**
     * Class providing core Lattice Boltzmann functionality.
     * Implements the IteratedAction interface.
     */
    template<class TRAITS = hemelb::Traits<>>
    class LBM : public LBMBase
    {
    public:
        //! Instantiation type
        using Traits = TRAITS;
    private:
        using LatticeType = typename Traits::Lattice;
        // Use the kernel specified through the build system. This will select one of the above classes.
        using LBKernel = typename Traits::Kernel;

        using tMidFluidCollision = typename Traits::Streamer;
        // Use the wall boundary condition specified through the build system.
        using tWallCollision = typename Traits::WallBoundary;
        // Use the inlet BC specified by the build system
        using tInletCollision = typename Traits::InletBoundary;
        // Use the outlet BC specified by the build system
        using tOutletCollision = typename Traits::OutletBoundary;
        // And again but for sites that are both in-/outlet and wall
        using tInletWallCollision = typename Traits::WallInletBoundary;
        using tOutletWallCollision = typename Traits::WallOutletBoundary;

    public:
        using LBMBase::LBMBase;
        ~LBM() override = default;


        void PreSend() override; ///< part of IteratedAction interface.
        void PreReceive() override; ///< part of IteratedAction interface.
        void PostReceive() override; ///< part of IteratedAction interface.


        void SetInitialConditions(lb::InitialCondition const& ic_conf, const net::IOCommunicator& ioComms) override;

    private:

        void InitCollisions() override;

        /// Collision objects
        /// We hold these in a tuple to allow heterogeneous iteration over them.
        // TODO: these should probably be optional but we can't because some of the collisions aren't copyable
        using CollisionTuple = std::tuple<
                std::unique_ptr<tMidFluidCollision>,
                std::unique_ptr<tWallCollision>,
                std::unique_ptr<tInletCollision>,
                std::unique_ptr<tOutletCollision>,
                std::unique_ptr<tInletWallCollision>,
                std::unique_ptr<tOutletWallCollision>
        >;
        CollisionTuple mCollisions;

        void StreamAndCollide(streamer auto& s, const site_t beginIndex, const site_t endIndex)
        {
            s.StreamAndCollide(beginIndex, endIndex, &mParams, *mLatDat, propertyCache);
        }

        void PostStep(streamer auto& s, const site_t beginIndex, const site_t endIndex)
        {
            s.PostStep(beginIndex, endIndex, &mParams, *mLatDat, propertyCache);
        }

    };
}
#endif // HEMELB_LB_LB_H
