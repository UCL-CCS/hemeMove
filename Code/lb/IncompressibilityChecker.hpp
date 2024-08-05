// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_INCOMPRESSIBILITYCHECKER_HPP
#define HEMELB_LB_INCOMPRESSIBILITYCHECKER_HPP

#include "lb/IncompressibilityChecker.h"

#include "hassert.h"
#include "reporting/Timers.h"

namespace hemelb::lb
{
    template<class BroadcastPolicy>
    IncompressibilityChecker<BroadcastPolicy>::DensityTracker::DensityTracker() :
        allocatedHere(true)
    {
      densitiesArray = new distribn_t[DENSITY_TRACKER_SIZE];

      //! @todo #23 do we have a policy on floating point constants?
      densitiesArray[MIN_DENSITY] = DBL_MAX;
      densitiesArray[MAX_DENSITY] = -DBL_MAX;
      densitiesArray[MAX_VELOCITY_MAGNITUDE] = 0.0;
    }

    template<class BroadcastPolicy>
    IncompressibilityChecker<BroadcastPolicy>::DensityTracker::DensityTracker(
        distribn_t* const densityValues) :
        densitiesArray(densityValues), allocatedHere(false)
    {
    }

    template<class BroadcastPolicy>
    IncompressibilityChecker<BroadcastPolicy>::DensityTracker::~DensityTracker()
    {
      if (allocatedHere)
      {
        delete[] densitiesArray;
      }
    }

    template<class BroadcastPolicy>
    void IncompressibilityChecker<BroadcastPolicy>::DensityTracker::operator=(
        const DensityTracker& newValues)
    {
      for (unsigned trackerEntry = 0; trackerEntry < DENSITY_TRACKER_SIZE; trackerEntry++)
      {
        densitiesArray[trackerEntry] = newValues.GetDensitiesArray()[trackerEntry];
      }
    }

    template<class BroadcastPolicy>
    distribn_t& IncompressibilityChecker<BroadcastPolicy>::DensityTracker::operator[](
        DensityTrackerIndices densityIndex) const
    {
      return densitiesArray[densityIndex];
    }

    template<class BroadcastPolicy>
    distribn_t* IncompressibilityChecker<BroadcastPolicy>::DensityTracker::GetDensitiesArray() const
    {
      return densitiesArray;
    }

    template<class BroadcastPolicy>
    void IncompressibilityChecker<BroadcastPolicy>::DensityTracker::UpdateDensityTracker(
        const DensityTracker& newValues)
    {
      if (newValues[MIN_DENSITY] < densitiesArray[MIN_DENSITY])
      {
        densitiesArray[MIN_DENSITY] = newValues[MIN_DENSITY];
      }
      if (newValues[MAX_DENSITY] > densitiesArray[MAX_DENSITY])
      {
        densitiesArray[MAX_DENSITY] = newValues[MAX_DENSITY];
      }
      if (newValues[MAX_VELOCITY_MAGNITUDE] > densitiesArray[MAX_VELOCITY_MAGNITUDE])
      {
        densitiesArray[MAX_VELOCITY_MAGNITUDE] = newValues[MAX_VELOCITY_MAGNITUDE];
      }
    }

    template<class BroadcastPolicy>
    void IncompressibilityChecker<BroadcastPolicy>::DensityTracker::UpdateDensityTracker(
        distribn_t newDensity, distribn_t newVelocityMagnitude)
    {
      if (newDensity < densitiesArray[MIN_DENSITY])
      {
        densitiesArray[MIN_DENSITY] = newDensity;
      }
      if (newDensity > densitiesArray[MAX_DENSITY])
      {
        densitiesArray[MAX_DENSITY] = newDensity;
      }
      if (newVelocityMagnitude > densitiesArray[MAX_VELOCITY_MAGNITUDE])
      {
        densitiesArray[MAX_VELOCITY_MAGNITUDE] = newVelocityMagnitude;
      }
    }

    template<class BroadcastPolicy>
    IncompressibilityChecker<BroadcastPolicy>::IncompressibilityChecker(
            const geometry::Domain * latticeData, net::Net* net, SimulationState* simState,
            lb::MacroscopicPropertyCache& propertyCache, reporting::Timers& timings,
            distribn_t maximumRelativeDensityDifferenceAllowed) :
        BroadcastPolicy(net, simState, SPREADFACTOR), mLatDat(latticeData),
            propertyCache(propertyCache), mSimState(simState), timings(timings),
            maximumRelativeDensityDifferenceAllowed(maximumRelativeDensityDifferenceAllowed),
            globalDensityTracker(nullptr)
    {
      /*
       *  childrenDensitiesSerialised must be initialised to something sensible since ReceiveFromChildren won't
       *  fill it in completely unless the logarithm base SPREADFACTOR of the number of processes is an integer.
       */
      for (unsigned leaf_index = 0; leaf_index < SPREADFACTOR; leaf_index++)
      {
        unsigned offset = leaf_index * DensityTracker::DENSITY_TRACKER_SIZE;
        for (unsigned tracker_index = 0; tracker_index < DensityTracker::DENSITY_TRACKER_SIZE;
            tracker_index++)
        {
          switch (tracker_index)
          {
            case DensityTracker::MIN_DENSITY:
            case DensityTracker::MAX_DENSITY:
              childrenDensitiesSerialised[offset + tracker_index] = REFERENCE_DENSITY;
              break;
            case DensityTracker::MAX_VELOCITY_MAGNITUDE:
              childrenDensitiesSerialised[offset + tracker_index] = 0.0;
              break;
            default:
              // This should never trip. It only occurs when a new entry is added to the density tracker and
              // no suitable initialisation is provided.
              HASSERT(false);
          }
        }
      }

    }

    template<class BroadcastPolicy>
    distribn_t IncompressibilityChecker<BroadcastPolicy>::GetGlobalSmallestDensity() const
    {
      HASSERT(AreDensitiesAvailable());
      return (*globalDensityTracker)[DensityTracker::MIN_DENSITY];
    }

    template<class BroadcastPolicy>
    distribn_t IncompressibilityChecker<BroadcastPolicy>::GetGlobalLargestDensity() const
    {
      HASSERT(AreDensitiesAvailable());
      return (*globalDensityTracker)[DensityTracker::MAX_DENSITY];
    }

    template<class BroadcastPolicy>
    double IncompressibilityChecker<BroadcastPolicy>::GetMaxRelativeDensityDifference() const
    {
      distribn_t maxDensityDiff = GetGlobalLargestDensity() - GetGlobalSmallestDensity();
      HASSERT(maxDensityDiff >= 0.0);
      return maxDensityDiff / REFERENCE_DENSITY;
    }

    template<class BroadcastPolicy>
    double IncompressibilityChecker<BroadcastPolicy>::GetMaxRelativeDensityDifferenceAllowed() const
    {
      return maximumRelativeDensityDifferenceAllowed;
    }

    template<class BroadcastPolicy>
    void IncompressibilityChecker<BroadcastPolicy>::PostReceiveFromChildren(
        unsigned long splayNumber)
    {
      timings.monitoring().Start();

      for (int childIndex = 0; childIndex < (int) SPREADFACTOR; childIndex++)
      {
        DensityTracker childDensities(&childrenDensitiesSerialised[childIndex
            * DensityTracker::DENSITY_TRACKER_SIZE]);

        upwardsDensityTracker.UpdateDensityTracker(childDensities);
      }

      timings.monitoring().Stop();
    }

    template<class BroadcastPolicy>
    void IncompressibilityChecker<BroadcastPolicy>::PostSendToParent(unsigned long splayNumber)
    {
      timings.monitoring().Start();

      for (site_t i = 0; i < mLatDat->GetLocalFluidSiteCount(); i++)
      {
        upwardsDensityTracker.UpdateDensityTracker(propertyCache.densityCache.Get(i),
                                                   propertyCache.velocityCache.Get(i).GetMagnitude());
      }

      timings.monitoring().Stop();
    }

    template<class BroadcastPolicy>
    void IncompressibilityChecker<BroadcastPolicy>::ProgressFromChildren(unsigned long splayNumber)
    {
      this->ReceiveFromChildren(childrenDensitiesSerialised, DensityTracker::DENSITY_TRACKER_SIZE);
    }

    template<class BroadcastPolicy>
    void IncompressibilityChecker<BroadcastPolicy>::ProgressFromParent(unsigned long splayNumber)
    {
      this->ReceiveFromParent(downwardsDensityTracker.GetDensitiesArray(),
                              DensityTracker::DENSITY_TRACKER_SIZE);
    }

    template<class BroadcastPolicy>
    void IncompressibilityChecker<BroadcastPolicy>::ProgressToChildren(unsigned long splayNumber)
    {
      this->SendToChildren(downwardsDensityTracker.GetDensitiesArray(),
                           DensityTracker::DENSITY_TRACKER_SIZE);
    }

    template<class BroadcastPolicy>
    void IncompressibilityChecker<BroadcastPolicy>::ProgressToParent(unsigned long splayNumber)
    {
      this->SendToParent(upwardsDensityTracker.GetDensitiesArray(),
                         DensityTracker::DENSITY_TRACKER_SIZE);
    }

    template<class BroadcastPolicy>
    void IncompressibilityChecker<BroadcastPolicy>::TopNodeAction()
    {
      downwardsDensityTracker = upwardsDensityTracker;
    }

    template<class BroadcastPolicy>
    void IncompressibilityChecker<BroadcastPolicy>::Effect()
    {
      globalDensityTracker = &downwardsDensityTracker;
    }

    template<class BroadcastPolicy>
    bool IncompressibilityChecker<BroadcastPolicy>::AreDensitiesAvailable() const
    {
      return (globalDensityTracker != nullptr);
    }

    template<class BroadcastPolicy>
    bool IncompressibilityChecker<BroadcastPolicy>::IsDensityDiffWithinRange() const
    {
      return (GetMaxRelativeDensityDifference() < maximumRelativeDensityDifferenceAllowed);
    }

    template<class BroadcastPolicy>
    void IncompressibilityChecker<BroadcastPolicy>::Report(reporting::Dict& dictionary)
    {
      if (AreDensitiesAvailable() && !IsDensityDiffWithinRange())
      {
        reporting::Dict incomp = dictionary.AddSectionDictionary("DENSITIES");
        incomp.SetFormattedValue("ALLOWED", "%.1f%%", GetMaxRelativeDensityDifferenceAllowed() * 100);
        incomp.SetFormattedValue("ACTUAL", "%.1f%%", GetMaxRelativeDensityDifference() * 100);
      }
    }

    template<class BroadcastPolicy>
    double IncompressibilityChecker<BroadcastPolicy>::GetGlobalLargestVelocityMagnitude() const
    {
      HASSERT(AreDensitiesAvailable());
      return (*globalDensityTracker)[DensityTracker::MAX_VELOCITY_MAGNITUDE];
    }
}

#endif /* HEMELB_LB_INCOMPRESSIBILITYCHECKER_HPP */
