// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_LBMPARAMETERS_H
#define HEMELB_LB_LBMPARAMETERS_H

#include <cmath>
#include <vector>

#include "constants.h"

namespace hemelb::geometry {
    class Domain;
    namespace neighbouring {
        class NeighbouringDataManager;
    }
}

namespace hemelb::lb
{
    class BoundaryValues;

    class LbmParameters
    {
        inline void CalcDerivedParams() {
            tau = 0.5 + (timeStep * eta / fluidDensity) / (Cs2 * voxelSize * voxelSize);
            omega = -1.0 / tau;
            stressParameter = (1.0 - 1.0 / (2.0 * tau)) / std::sqrt(2.0);
            beta = -1.0 / (2.0 * tau);
        }
      public:
        inline LbmParameters()
        {
            CalcDerivedParams();
        }

        inline LbmParameters(PhysicalTime timeStepSeconds,
                             PhysicalDistance voxelSizeMetres) :
                timeStep(timeStepSeconds),
                voxelSize(voxelSizeMetres)
        {
            CalcDerivedParams();
        }
        inline LbmParameters(PhysicalTime timeStepSeconds,
                             PhysicalDistance voxelSizeMetres,
                             PhysicalDensity fluidDensityKgm3
        ) :
                timeStep(timeStepSeconds),
                voxelSize(voxelSizeMetres),
                fluidDensity(fluidDensityKgm3)
        {
            CalcDerivedParams();
        }
        inline LbmParameters(PhysicalTime timeStepSeconds,
                             PhysicalDistance voxelSizeMetres,
                             PhysicalDensity fluidDensityKgm3,
                             PhysicalDynamicViscosity newtonianViscosityPas) :
                timeStep(timeStepSeconds),
                voxelSize(voxelSizeMetres),
                fluidDensity(fluidDensityKgm3),
                eta(newtonianViscosityPas)
        {
            CalcDerivedParams();
        }
        inline const PhysicalTime& GetTimeStep() const
        {
          return timeStep;
        }

        inline const PhysicalDistance& GetVoxelSize() const
        {
          return voxelSize;
        }

        inline const PhysicalDensity& GetFluidDensity() const
        {
          return fluidDensity;
        }

        inline const PhysicalDynamicViscosity& GetEta() const
        {
          return eta;
        }

        inline const distribn_t& GetOmega() const
        {
          return omega;
        }

        inline const distribn_t& GetTau() const
        {
          return tau;
        }

        inline const distribn_t& GetStressParameter() const
        {
          return stressParameter;
        }

        inline const distribn_t& GetBeta() const
        {
          return beta;
        }

      private:
        PhysicalTime timeStep = 1; // seconds
        PhysicalDistance voxelSize = 1; // metres
        PhysicalDensity fluidDensity = DEFAULT_FLUID_DENSITY_Kg_per_m3; // kg m^-3
        PhysicalDynamicViscosity eta = DEFAULT_FLUID_VISCOSITY_Pas; // Newtonian viscosity of fluid in Pa s
        distribn_t tau;
        distribn_t omega;
        distribn_t stressParameter;
        distribn_t beta; ///< Viscous dissipation in ELBM
    };
/**
     * InitParams: struct for passing variables into streaming, collision and kernel operators
     * to initialise them.
     *
     * When a newly-developed kernel, collider or streamer requires extra parameters to be
     * passed in for initialisation, it's annoying to have to change the constructors in
     * multiple places to make them all consistent (so that higher-up code can seamlessly
     * construct one kind or another).
     *
     * Instead, new parameters can be added to this single object, which should be the only
     * constructor argument used by any kernel / collision / streaming implementation.
     */
    struct InitParams {
    public:

        // Assume the first site to be used in the kernel is the first site in the core, unless otherwise specified
        InitParams() = default;

        // Each streamer is responsible for updating certain types of sites. These are arranged such they are largely
        // contiguous in memory (the local contiguous site id). This data structure refers to which of those are handled
        // by the current streamer. These are given as a collection of contiguous site ids, running from e.g.
        // siteRanges[0].first to siteRanges[0].second-1 (inclusive).
        std::vector<std::pair<site_t, site_t> > siteRanges;

        // The array with the imposed density at each boundary.
        BoundaryValues *boundaryObject;

        // The lattice data object. Currently only used for accessing the boundary id
        // of each site next to an inlet or an outlet.
        const geometry::Domain *latDat;

        // The LB parameters object. Currently only used in LBGKNN to access the current
        // time step.
        const LbmParameters *lbmParams;

        // The neighbouring data manager, for kernels / collisions / streamers that
        // require data from other cores.
        geometry::neighbouring::NeighbouringDataManager *neighbouringDataManager;
    };
}

#endif //HEMELB_LB_LBMPARAMETERS_H
