// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_STABILITYTESTER_H
#define HEMELB_LB_STABILITYTESTER_H

#include "net/PhasedBroadcastRegular.h"
#include "geometry/Domain.h"
#include "geometry/FieldData.h"
#include "configuration/MonitoringConfig.h"
#include "reporting/Timers.h"

namespace hemelb::lb
{
    /**
     * Class to repeatedly assess the stability of the simulation, using the PhasedBroadcast
     * interface.
     *
     * PhasedBroadcastRegular is used because we know in advance which iterations we will
     * want to communicate on. The default parameters suffice: no initial action is necessary
     * because we can assess the stability just before communicating (it doesn't have to happen
     * at the same time on all nodes), only one communication is needed between depths, which
     * can't overlap. We go down the tree to pass the overall stability to all nodes, and we go up
     * the tree to compose the local stability for all nodes to discover whether the simulation as
     * a whole is stable.
     */
    class StabilityTester: public net::PhasedBroadcastRegular<> {
    protected:
        /**
         * Slightly arbitrary spread factor for the tree.
         */
        static constexpr unsigned SPREADFACTOR = 10;

        std::shared_ptr<const geometry::FieldData> mLatDat;

        /**
         * Stability value of this node and its children to propagate upwards.
         */
        int mUpwardsStability;
        /**
         * Stability value as understood by the root node, to pass downwards.
         */
        int mDownwardsStability;
        /**
         * Array for storing the passed-up stability values from child nodes.
         */
        std::array<int, SPREADFACTOR> mChildrensStability;
        /**
         * Pointer to the simulation state used in the rest of the simulation.
         */
        lb::SimulationState* mSimState;

        /** Timing object. */
        reporting::Timers& timings;

        /** Object containing the user-provided configuration for this class */
        hemelb::configuration::MonitoringConfig testerConfig;

    public:
        StabilityTester(std::shared_ptr<const geometry::FieldData> iLatDat, net::Net* net,
                        SimulationState* simState, reporting::Timers& timings,
                        const hemelb::configuration::MonitoringConfig& testerConfig);

        bool ShouldTerminateWhenConverged() const;
        /**
         * Override the reset method in the base class, to reset the stability variables.
         */
        void Reset();
    protected:
        /**
         * Override the methods from the base class to propagate data from the root, and
         * to send data about this node and its childrens' stabilities up towards the root.
         */
        void ProgressFromChildren(unsigned long splayNumber) override;

        void ProgressFromParent(unsigned long splayNumber) override;

        void ProgressToChildren(unsigned long splayNumber) override;

        void ProgressToParent(unsigned long splayNumber) override;

        // PostSendToParent in Impl below

        /// Take the combined stability information (an int, with a value of hemelb::lb::Unstable
        /// if any child node is unstable) and start passing it back down the tree.
        void TopNodeAction() override;
        /// Override the method from the base class to use the data from child nodes.
        void PostReceiveFromChildren(unsigned long splayNumber) override;

        /**
         * Apply the stability value sent by the root node to the simulation logic.
         */
        void Effect() override;
    };

    template<class LatticeType>
    class StabilityTesterImpl : public StabilityTester
    {
    public:
        using StabilityTester::StabilityTester;

    protected:
        /**
         * The algorithm that checks distribution function convergence must be run in this
         * method rather than in ProgressToParent to make sure that the current timestep has
         * finished streaming.
         *
         * @param splayNumber
         */
        void PostSendToParent(unsigned long splayNumber) override
        {
          timings.monitoring().Start();
          // No need to bother testing out local lattice points if we're going to be
          // sending up a 'Unstable' value anyway.
          if (mUpwardsStability != Unstable)
          {
            bool unconvergedSitePresent = false;

            for (site_t i = 0; i < mLatDat->GetDomain().GetLocalFluidSiteCount(); i++)
            {
              for (unsigned int l = 0; l < LatticeType::NUMVECTORS; l++)
              {
                distribn_t value = *mLatDat->GetFNew(i * LatticeType::NUMVECTORS + l);

                // Note that by testing for value > 0.0, we also catch stray NaNs.
                if (! (value > 0.0))
                {
                  mUpwardsStability = Unstable;
                  break;
                }
              }
              ///@todo: If we refactor the previous loop out, we can get away with a single break statement
              if (mUpwardsStability == Unstable)
              {
                break;
              }

              if (testerConfig.doConvergenceCheck)
              {
                distribn_t relativeDifference =
                    ComputeRelativeDifference(mLatDat->GetFNew<LatticeType>(i),
                                              mLatDat->GetSite(i).GetFOld<LatticeType>());

                if (relativeDifference > testerConfig.convergenceRelativeTolerance)
                {
                  // The simulation is stable but hasn't converged in the whole domain yet.
                  unconvergedSitePresent = true;
                }
              }
            }

            switch (mUpwardsStability)
            {
              case UndefinedStability:
              case Stable:
              case StableAndConverged:
                mUpwardsStability = (testerConfig.doConvergenceCheck && !unconvergedSitePresent) ?
                  StableAndConverged :
                  Stable;
                break;
              case Unstable:
                break;
            }
          }

          timings.monitoring().Stop();
        }

        using span = typename LatticeType::const_span;
        /**
         * Computes the relative difference between the densities at the beginning and end of a
         * timestep, i.e. |(rho_new - rho_old) / (rho_old - rho_0)|.
         *
         * @param fNew Distribution function after stream and collide, i.e. solution of the current timestep.
         * @param fOld Distribution function at the end of the previous timestep.
         * @return relative difference between the densities computed from fNew and fOld.
         */
        inline double ComputeRelativeDifference(span fNew,
                                                span fOld) const
        {
            distribn_t newDensity;
            LatticeMomentum newMomentum;
            LatticeType::CalculateDensityAndMomentum(fNew,
                                                     newDensity,
                                                     newMomentum);

            distribn_t oldDensity;
            LatticeMomentum oldMomentum;
            LatticeType::CalculateDensityAndMomentum(fOld,
                                                     oldDensity,
                                                     oldMomentum);

            if (std::holds_alternative<extraction::source::Velocity>(testerConfig.convergenceVariable))
            {
                auto diff_vel = newMomentum / newDensity - oldMomentum / oldDensity;
                auto absoluteError = diff_vel.GetMagnitude();
                return absoluteError / testerConfig.convergenceReferenceValue;
            } else {
                throw Exception() << "Convergence check based on requested variable currently not available";
            }
        }
    };
}

#endif /* HEMELB_LB_STABILITYTESTER_H */
