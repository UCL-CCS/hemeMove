// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_COLLISIONS_NONZEROVELOCITYEQUILIBRIUMFIXEDDENSITY_H
#define HEMELB_LB_COLLISIONS_NONZEROVELOCITYEQUILIBRIUMFIXEDDENSITY_H

#include "lb/concepts.h"
#include "lb/iolets/BoundaryValues.h"
#include "geometry/Site.h"

namespace hemelb::lb
{

    template<kernel_type K>
    class NonZeroVelocityEquilibriumFixedDensity
    {
    public:
        using KernelType = K;
        using LatticeType = typename KernelType::LatticeType;
        using VarsType = typename KernelType::VarsType;

        NonZeroVelocityEquilibriumFixedDensity(InitParams& initParams) :
                kernel(initParams), boundaryObject(initParams.boundaryObject)
        {
        }

        void CalculatePreCollision(VarsType& hydroVars,
                                   const geometry::Site<geometry::Domain>& site)
        {
            LatticeType::CalculateDensityAndMomentum(hydroVars.f,
                                                     hydroVars.density,
                                                     hydroVars.momentum);

            // Externally impose a density. Keep a record of the old one so we can scale the
            // momentum vector.
            distribn_t previousDensity = hydroVars.density;
            hydroVars.density = boundaryObject->GetBoundaryDensity(site.GetIoletId());

            hydroVars.momentum *= (hydroVars.density / previousDensity);

            kernel.CalculateFeq(hydroVars, site.GetIndex());
        }

        void Collide(const LbmParameters* lbmParams,
                     VarsType& iHydroVars)
        {
            for (Direction direction = 0; direction < LatticeType::NUMVECTORS; ++direction)
            {
                iHydroVars.SetFPostCollision(direction, iHydroVars.GetFEq()[direction]);
            }
        }

    private:
        KernelType kernel;
        BoundaryValues* boundaryObject;
    };
}

#endif
