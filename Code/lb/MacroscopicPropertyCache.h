// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_MACROSCOPICPROPERTYCACHE_H
#define HEMELB_LB_MACROSCOPICPROPERTYCACHE_H

#include <vector>
#include "lb/SimulationState.h"
#include "units.h"
#include "util/RefreshableCache.hpp"
#include "util/Matrix3D.h"

namespace hemelb
{
  namespace geometry
  {
    class Domain;
  }
  namespace lb
  {
    class MacroscopicPropertyCache
    {
      public:
        /**
         * Constructor, the only way to create this object.
         * @param simState The simulation state, so that the cache knows when it is out of date.
         * @param fluidSiteCount Number of fluid sites, so the cache knows how big it needs to be.
         * @return
         */
        MacroscopicPropertyCache(const SimulationState& simState,
                                 const geometry::Domain& latticeData);

        /**
         * Reset the list of cache types required to be none of them.
         */
        void ResetRequirements();

        /**
         * Returns the number of sites cached.
         * @return
         */
        site_t GetSiteCount() const;

        /**
         * The cache of densities for each fluid site on this core.
         */
        util::RefreshableCache<distribn_t> densityCache;

        /**
         * The cache of velocities for each fluid site on this core.
         */
        util::RefreshableCache<util::Vector3D<distribn_t> > velocityCache;

        /**
         * The cache of wall shear stress magnitudes for each fluid site on this core.
         */
        util::RefreshableCache<distribn_t> wallShearStressMagnitudeCache;

        /**
         * The cache of Von Mises stresses for each fluid site on this core.
         */
        util::RefreshableCache<distribn_t> vonMisesStressCache;

        /**
         * The cache of shear rates for each fluid site on this core.
         */
        util::RefreshableCache<distribn_t> shearRateCache;

        /**
         * The cache of stress tensors for each fluid site on this core.
         */
        util::RefreshableCache<util::Matrix3D> stressTensorCache;

        /**
         * The cache of traction vectors for each fluid site on this core.
         */
        util::RefreshableCache<util::Vector3D<LatticeStress> > tractionCache;

        /**
         * The cache of projections of the traction vectors on the tangential plane of each fluid site on this core.
         */
        util::RefreshableCache<util::Vector3D<LatticeStress> > tangentialProjectionTractionCache;

        /**
         * The cache of the velocity distributions for each fluid site on this core.
         */
        util::RefreshableCache<util::Vector3D<LatticeStress> > velDistributionsCache;

      private:
        /**
         * The state of the simulation, including the number of timesteps passed.
         */
        // const SimulationState& simulationState;
        /**
         * The number of sites.
         */
        site_t siteCount;
    };
  }
}

#endif /* HEMELB_LB_MACROSCOPICPROPERTYCACHE_H */
