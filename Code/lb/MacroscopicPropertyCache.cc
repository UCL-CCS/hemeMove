// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "lb/MacroscopicPropertyCache.h"

#include "geometry/Domain.h"

namespace hemelb::lb
{
    MacroscopicPropertyCache::MacroscopicPropertyCache(const SimulationState& simState,
                                                       const geometry::Domain& latticeData) :
      densityCache(simState, latticeData.GetLocalFluidSiteCount()),
      velocityCache(simState, latticeData.GetLocalFluidSiteCount()),
      wallShearStressMagnitudeCache(simState, latticeData.GetLocalFluidSiteCount()),
      vonMisesStressCache(simState, latticeData.GetLocalFluidSiteCount()),
      shearRateCache(simState, latticeData.GetLocalFluidSiteCount()),
      stressTensorCache(simState, latticeData.GetLocalFluidSiteCount()),
      tractionCache(simState, latticeData.GetLocalFluidSiteCount()),
      tangentialProjectionTractionCache(simState, latticeData.GetLocalFluidSiteCount()),
      velDistributionsCache(simState, latticeData.GetLocalFluidSiteCount()),
      siteCount(latticeData.GetLocalFluidSiteCount())
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
      velDistributionsCache.UnsetRefreshFlag();
    }

    site_t MacroscopicPropertyCache::GetSiteCount() const
    {
      return siteCount;
    }
}

