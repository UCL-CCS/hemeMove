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
                                                       const geometry::LatticeData& latticeData) :
      densityCache(simState, latticeData.GetLocalFluidSiteCount()),
          velocityCache(simState, latticeData.GetLocalFluidSiteCount()),
          shearStressCache(simState, latticeData.GetLocalFluidSiteCount()),
          vonMisesStressCache(simState, latticeData.GetLocalFluidSiteCount()),
          shearRateCache(simState, latticeData.GetLocalFluidSiteCount()),
          simulationState(simState),
          siteCount(latticeData.GetLocalFluidSiteCount())
    {
      std::cout << "MacroscopicPropertyCache created of size " << latticeData.GetLocalFluidSiteCount() << std::endl;
      ResetRequirements();
    }

    void MacroscopicPropertyCache::ResetRequirements()
    {
      densityCache.UnsetRefreshFlag();
      velocityCache.UnsetRefreshFlag();
      vonMisesStressCache.UnsetRefreshFlag();
      shearStressCache.UnsetRefreshFlag();
      shearRateCache.UnsetRefreshFlag();
    }

    site_t MacroscopicPropertyCache::GetSiteCount() const
    {
      return siteCount;
    }
  }
}

