// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_TESTS_LB_LBTESTSHELPER_H
#define HEMELB_TESTS_LB_LBTESTSHELPER_H

#include <cmath>
#include <string>

#include <catch2/catch.hpp>

#include "constants.h"
#include "lb/concepts.h"
#include "lb/HFunction.h"
#include "lb/MacroscopicPropertyCache.h"
#include "lb/lattices/D3Q15.h"

#include "tests/helpers/ApproxVector.h"
#include "tests/helpers/FourCubeLatticeData.h"

namespace hemelb::tests::LbTestsHelper
{
    template <lb::lattice_type LatticeType>
    auto ZeroArray() {
        std::array<distribn_t, LatticeType::NUMVECTORS> ans;
        std::fill(ans.begin(), ans.end(), 0.0);
        return ans;
    }

    template<lb::lattice_type Lattice>
    using const_span = typename Lattice::const_span;
    template<lb::lattice_type Lattice>
    using mut_span = typename Lattice::mut_span;

    using Vec = util::Vector3D<distribn_t>;

    template<lb::lattice_type Lattice>
    void CalculateRhoMomentum(const_span<Lattice> f, distribn_t& rho,  Vec& v)
    {
        rho = 0.0;
        v = {0.0, 0.0, 0.0};
        for (unsigned int ii = 0; ii < Lattice::NUMVECTORS; ++ii)
        {
            rho += f[ii];
            v += f[ii] * Lattice::VECTORS[ii];
        }
    }

    template<lb::lattice_type Lattice>
    Vec CalculateMomentum(const_span<Lattice> f)
    {
        auto v = Vec::Zero();
        for (unsigned int ii = 0; ii < Lattice::NUMVECTORS; ++ii) {
            v += f[ii] * Vec{Lattice::VECTORS[ii]};
        }
        return v;
    }

    template<lb::lattice_type Lattice>
    void CalculateAnsumaliEntropicEqmF(distribn_t density,
                                       const Vec& momentum,
                                       mut_span<Lattice> f)
    {
        // Calculate velocity.
        auto velocity = momentum / density;

        distribn_t B[3];
        distribn_t term1[3];
        distribn_t term2[3];
        for (int ii = 0; ii < 3; ++ii)
        {
            B[ii] = sqrt(1.0 + 3.0 * velocity[ii] * velocity[ii]);
            term1[ii] = 2.0 - B[ii];
            term2[ii] = (2.0 * velocity[ii] + B[ii]) / (1.0 - velocity[ii]);
        }

        for (unsigned int ii = 0; ii < Lattice::NUMVECTORS; ++ii)
        {
            f[ii] = density * Lattice::EQMWEIGHTS[ii];

            f[ii] *= term1[0] * pow(term2[0], Lattice::CXD[ii]);
            f[ii] *= term1[1] * pow(term2[1], Lattice::CYD[ii]);
            f[ii] *= term1[2] * pow(term2[2], Lattice::CZD[ii]);
        }
    }

    template<lb::lattice_type Lattice>
    void CalculateLBGKEqmF(distribn_t density,
                           const util::Vector3D<distribn_t>& momentum,
                           mut_span<Lattice> f)
    {
        for (unsigned int ii = 0; ii < Lattice::NUMVECTORS; ++ii)
        {
            // Calculate the dot-product of the velocity with the direction vector.
            distribn_t vSum = Dot(momentum, Lattice::VECTORS[ii].template as<double>());

            // Calculate the squared magnitude of the velocity.
            distribn_t momentum2Sum = momentum.GetMagnitudeSquared();

            // F eqm = density proportional component...
            f[ii] = density;

            // ... - v^2 component...
            f[ii] -= ( (3.0 / 2.0) * momentum2Sum / density);

            // ... + v^1 component
            f[ii] += 3.0 * vSum + (9.0 / 2.0) * vSum * vSum / density;

            // Multiply by eqm weight.
            f[ii] *= Lattice::EQMWEIGHTS[ii];
        }
    }

    template<lb::lattice_type Lattice>
    void CalculateEntropicCollision(const_span<Lattice> f,
                                    const_span<Lattice> f_eq,
                                    distribn_t tau,
                                    distribn_t beta,
                                    mut_span<Lattice> f_collided)
    {
        lb::HFunction<Lattice> HFunc(f, f_eq);

        distribn_t alpha = util::NewtonRaphson(HFunc, 2.0, 1.0E-100);

        for (unsigned int ii = 0; ii < Lattice::NUMVECTORS; ++ii)
        {
            f_collided[ii] = f[ii] + alpha * beta * (f[ii] - f_eq[ii]);
        }
    }

    template<lb::lattice_type Lattice>
    void CalculateLBGKCollision(const_span<Lattice> f,
                                const_span<Lattice> f_eq,
                                distribn_t omega,
                                mut_span<Lattice> f_collided)
    {
        for (unsigned int ii = 0; ii < Lattice::NUMVECTORS; ++ii)
        {
            f_collided[ii] = f[ii] + omega * (f[ii] - f_eq[ii]);
        }
    }

    //using DISTS = distribn_t[lb::D3Q15::NUMVECTORS];
    template<typename VarsType>
    void CompareHydros(const distribn_t expectedDensity,
                       const LatticeMomentum& expectedMomentum,
                       typename VarsType::const_span expectedFEq,
                       VarsType& hydroVars,
                       distribn_t allowedError)
    {
        // Compare density
        REQUIRE(Approx(expectedDensity).margin(allowedError) == hydroVars.density);
        // Compare momentum
        REQUIRE(ApproxVector<distribn_t>(expectedMomentum).Margin(allowedError) == hydroVars.momentum);

        // Compare equilibrium f
        for (unsigned int ii = 0; ii < lb::D3Q15::NUMVECTORS; ++ii) {
            REQUIRE(Approx(expectedFEq[ii]).margin(allowedError) == hydroVars.GetFEq()[ii]);
        }
    }

    template<lb::lattice_type LatticeType>
    void InitialiseAnisotropicTestData(site_t site, distribn_t* distribution)
    {
        for (unsigned int direction = 0; direction < LatticeType::NUMVECTORS; ++direction)
        {
            distribution[direction] = ((distribn_t) (direction + 1)) / 10.0 + ((distribn_t) (site)) / 100.0;
        }
    }

    template<lb::lattice_type LatticeType>
    void InitialiseAnisotropicTestData(FourCubeLatticeData& latticeData)
    {
        auto&& dom = latticeData.GetDomain();
        for (site_t site = 0; site < dom.GetLocalFluidSiteCount(); ++site)
        {
            distribn_t fOld[LatticeType::NUMVECTORS];
            InitialiseAnisotropicTestData<LatticeType> (site, (distribn_t*)fOld);
            latticeData.SetFOld<LatticeType>(site, fOld);
        }
    }


    /**
     * Updates a property cache for the macroscopic variables selected. This should have
     * identical behaviour to in UpdateSiteMinsAndMaxes in BaseStreamer where the cache
     * is normally populated.
     * @param latDat
     * @param cache
     * @param simState
     */
    template<lb::lattice_type Lattice>
    void UpdatePropertyCache(geometry::FieldData& latDat,
                             lb::MacroscopicPropertyCache& cache,
                             lb::SimulationState& simState)
    {
        const auto N = latDat.GetDomain().GetLocalFluidSiteCount();
        for (site_t site = 0; site < N; ++site)
        {
            distribn_t density, feq[Lattice::NUMVECTORS];
            util::Vector3D<distribn_t> momentum;
            util::Vector3D<distribn_t> velocity;

            Lattice::CalculateDensityMomentumFEq(latDat.GetSite(site).GetFOld<Lattice> (),
                                                 density,
                                                 momentum,
                                                 velocity,
                                                 feq);

            if (cache.densityCache.RequiresRefresh())
            {
                cache.densityCache.Put(site, density);
            }
            if (cache.velocityCache.RequiresRefresh())
            {
                cache.velocityCache.Put(site, velocity);
            }

            // TODO stress cache filling not yet implemented.
        }
    }

    /**
     * For every site in m_fieldData and every link crossing a boundary, set the distance to the
     * boundary to be wallDistance lattice length units
     *
     * @param latticeData Lattice object
     * @param wallDistance Distance to the wall in lattice units and in [0,1)
     */
    template<lb::lattice_type Lattice>
    void SetWallAndIoletDistances(FourCubeLatticeData& latticeData, distribn_t wallDistance)
    {
        assert(wallDistance >= 0);
        assert(wallDistance < 1);
        auto&& dom = latticeData.GetDomain();

        for (site_t siteIndex = 0; siteIndex < dom.GetTotalFluidSites(); siteIndex++)
        {
            for (Direction direction = 1; direction < Lattice::NUMVECTORS; direction++)
            {
                // -1 means that the a given link does not cross any boundary
                if (latticeData.GetSite(siteIndex).template GetWallDistance<Lattice> (direction) != (distribn_t) -1)
                {
                    dom.SetBoundaryDistance(siteIndex, direction, wallDistance);
                }
            }
        }
    }
}
#endif
