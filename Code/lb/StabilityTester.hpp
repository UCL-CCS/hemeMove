//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

#ifndef HEMELB_LB_STABILITYTESTER_HPP
#define HEMELB_LB_STABILITYTESTER_HPP

namespace hemelb
{
  namespace lb
  {
    /**
     * Class to repeatedly assess the stability of the simulation, using non-blocking collective
     */
    template<class LatticeType>
    StabilityTester<LatticeType>::StabilityTester(const geometry::LatticeData * iLatDat,
                                                  net::Net* net, SimulationState* simState,
                                                  reporting::Timers& timings,
                                                  bool checkForConvergence,
                                                  double relativeTolerance) :
        mLatDat(iLatDat), mSimState(simState), timings(timings),
            checkForConvergence(checkForConvergence), relativeTolerance(relativeTolerance),
            collectiveComm(net->GetCommunicator().Duplicate())
    {
      Reset();
    }

    /**
     * Override the reset method in the base class, to reset the stability variables.
     */
    template<class LatticeType>
    void StabilityTester<LatticeType>::Reset()
    {
      localStability = UndefinedStability;
      globalStability = UndefinedStability;
      mSimState->SetStability(UndefinedStability);
    }

    template<class LatticeType>
    bool StabilityTester<LatticeType>::CallAction(int action)
    {
      switch (static_cast<net::phased::steps::Step>(action))
      {
        case net::phased::steps::PreSend:
          PreSend();
          return true;
        case net::phased::steps::Send:
          Send();
          return true;
        case net::phased::steps::Wait:
          Wait();
          return true;
        case net::phased::steps::EndPhase:
          Effect();
          return true;
        default:
          return false;
      }
    }

    /**
     * Compute the local stability/convergence state.
     */
    template<class LatticeType>
    void StabilityTester<LatticeType>::PreSend(void)
    {
      timings[hemelb::reporting::Timers::monitoring].Start();
      bool unconvergedSitePresent = false;
      bool checkConvThisTimeStep = checkForConvergence;
      localStability = Stable;

      for (site_t i = 0; i < mLatDat->GetLocalFluidSiteCount(); i++)
      {
        for (unsigned int l = 0; l < LatticeType::NUMVECTORS; l++)
        {
          distribn_t value = *mLatDat->GetFNew(i * LatticeType::NUMVECTORS + l);

          // Note that by testing for value > 0.0, we also catch stray NaNs.
          if (! (value > 0.0))
          {
            localStability = Unstable;
            break;
          }
        }
        ///@todo: If we refactor the previous loop out, we can get away with a single break statement
        if (localStability == Unstable)
          break;

        if (checkConvThisTimeStep)
        {
          distribn_t relativeDifference = ComputeRelativeDifference(mLatDat->GetFNew(i
                                                                        * LatticeType::NUMVECTORS),
                                                                    mLatDat->GetSite(i).GetFOld<
                                                                        LatticeType>());

          if (relativeDifference > relativeTolerance)
          {
            // The simulation is stable but hasn't converged in the whole domain yet.
            unconvergedSitePresent = true;
            // Skip further tests for this timestep
            checkConvThisTimeStep = false;
          }
        }
      }

      if (localStability == Stable)
      {
        if (checkForConvergence && !unconvergedSitePresent)
        {
          localStability = StableAndConverged;
        }
      }
      timings[hemelb::reporting::Timers::monitoring].Stop();
    }

    /**
     * Initiate the collective.
     */
    template<class LatticeType>
    void StabilityTester<LatticeType>::Send(void)
    {
      // Begin collective.
      HEMELB_MPI_CALL(MPI_Iallreduce,
                      (static_cast<void*>(&localStability), static_cast<void*>(&globalStability), 1, net::MpiDataType(localStability), MPI_MIN, collectiveComm, &collectiveReq));
    }

    /**
     * Wait on the collectives to finish.
     */
    template<class LatticeType>
    void StabilityTester<LatticeType>::Wait(void)
    {
      timings[hemelb::reporting::Timers::mpiWait].Start();
      HEMELB_MPI_CALL(MPI_Wait, (&collectiveReq, MPI_STATUS_IGNORE));
      timings[hemelb::reporting::Timers::mpiWait].Stop();
    }

    /**
     * Computes the relative difference between the densities at the beginning and end of a
     * timestep, i.e. |(rho_new - rho_old) / (rho_old - rho_0)|.
     *
     * @param fNew Distribution function after stream and collide, i.e. solution of the current timestep.
     * @param fOld Distribution function at the end of the previous timestep.
     * @return relative difference between the densities computed from fNew and fOld.
     */
    template<class LatticeType>
    double StabilityTester<LatticeType>::ComputeRelativeDifference(const distribn_t* fNew,
                                                                   const distribn_t* fOld) const
    {
      distribn_t new_density = 0.;
      distribn_t old_density = 0.;
      for (unsigned int l = 0; l < LatticeType::NUMVECTORS; l++)
      {
        new_density += fNew[l];
        old_density += fOld[l];
      }

      // This is equivalent to REFERENCE_PRESSURE_mmHg in lattice units
      distribn_t ref_density = 1.0;

      if (old_density == ref_density)
      {
        // We want to avoid returning inf if the site is at pressure = REFERENCE_PRESSURE_mmHg
        return 0.0;
      }

      return fabs( (new_density - old_density) / (old_density - ref_density));
    }

    /**
     * Apply the stability value sent by the root node to the simulation logic.
     */
    template<class LatticeType>
    void StabilityTester<LatticeType>::Effect()
    {
      mSimState->SetStability((Stability) globalStability);
    }
  }
}

#endif /* HEMELB_LB_STABILITYTESTER_HPP */
