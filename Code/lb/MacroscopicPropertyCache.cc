// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#include "lb/MacroscopicPropertyCache.h"

namespace hemelb
{
  namespace lb
  {
    MacroscopicPropertyCache::MacroscopicPropertyCache(const SimulationState& simState,
                                                       site_t fluidSiteCount) :
      densityCache(simState, fluidSiteCount),
      velocityCache(simState, fluidSiteCount),
      wallShearStressMagnitudeCache(simState, fluidSiteCount),
      vonMisesStressCache(simState, fluidSiteCount),
      shearRateCache(simState, fluidSiteCount),
      stressTensorCache(simState, fluidSiteCount),
      tractionCache(simState, fluidSiteCount),
      tangentialProjectionTractionCache(simState, fluidSiteCount),
      simulationState(simState),
      siteCount(fluidSiteCount)
    {
      ResetRequirements();
    }

    void MacroscopicPropertyCache::ResetRequirements()
    {
      densityCache.UnsetRefreshFlag();
      velocityCache.UnsetRefreshFlag();
      vonMisesStressCache.UnsetRefreshFlag();
      wallShearStressMagnitudeCache.UnsetRefreshFlag();
      shearRateCache.UnsetRefreshFlag();
      stressTensorCache.UnsetRefreshFlag();
      tractionCache.UnsetRefreshFlag();
      tangentialProjectionTractionCache.UnsetRefreshFlag();
    }

    site_t MacroscopicPropertyCache::GetSiteCount() const
    {
      return siteCount;
    }
  }
}

