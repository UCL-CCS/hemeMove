// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_KERNELS_TRT_H
#define HEMELB_LB_KERNELS_TRT_H

#include <cstdlib>
#include "lb/HFunction.h"
#include "util/numerical.h"

namespace hemelb::lb
{
    /**
     * TRT: This class implements a two-relaxation time kernel.
     */
    template<lattice_type L>
    class TRT
    {
    public:
        using LatticeType = L;
        using VarsType = HydroVars<TRT>;

    private:
        static constexpr Direction FindZeroIndex() {
            for (Direction i = 0; i < LatticeType::NUMVECTORS; ++i) {
                Direction iBar = LatticeType::INVERSEDIRECTIONS[i];
                if (i == iBar)
                    return i;
            }
            return LatticeType::NUMVECTORS;
        }

        static constexpr Direction iZero = FindZeroIndex();
        static constexpr bool HasZero = iZero < LatticeType::NUMVECTORS;
        // Odd number => has a zero, so rounding down correct.
        static constexpr std::size_t NPAIRS = LatticeType::NUMVECTORS / 2;
        // Store the directions as pairs of opposites
        // (Note that the zero vector is its own opposite)
        using Opposites = std::pair<Direction, Direction>;
        static constexpr auto MakeOpposites() {
            std::array<Opposites, NPAIRS> ans;
            std::size_t j = 0;
            for (Direction i = 0; i < LatticeType::NUMVECTORS; ++i)
            {
                Direction iBar = LatticeType::INVERSEDIRECTIONS[i];
                if (iBar >= i) {
                    ans[j] = {i, iBar};
                    ++j;
                }
            }
            return ans;
        }
        static constexpr auto directionPairs = MakeOpposites();

    public:
        TRT(InitParams& initParams)
        {
        }

        void CalculateDensityMomentumFeq(VarsType& hydroVars, site_t index)
        {
            LatticeType::CalculateDensityMomentumFEq(hydroVars.f,
                                                     hydroVars.density,
                                                     hydroVars.momentum.x,
                                                     hydroVars.momentum.y,
                                                     hydroVars.momentum.z,
                                                     hydroVars.velocity.x,
                                                     hydroVars.velocity.y,
                                                     hydroVars.velocity.z,
                                                     hydroVars.f_eq.f);

            for (unsigned int ii = 0; ii < LatticeType::NUMVECTORS; ++ii)
            {
                hydroVars.f_neq.f[ii] = hydroVars.f[ii] - hydroVars.f_eq.f[ii];
            }
        }

        void CalculateFeq(VarsType& hydroVars, site_t index)
        {
            LatticeType::CalculateFeq(hydroVars.density,
                                      hydroVars.momentum.x,
                                      hydroVars.momentum.y,
                                      hydroVars.momentum.z,
                                      hydroVars.f_eq.f);

            for (unsigned int ii = 0; ii < LatticeType::NUMVECTORS; ++ii)
            {
                hydroVars.f_neq.f[ii] = hydroVars.f[ii] - hydroVars.f_eq.f[ii];
            }
        }

        void Collide(const LbmParameters* const lbmParams, VarsType& hydroVars)
        {
            // Note HemeLB defines omega = -1/ tau
            // Magic number determines the other relaxation time
            // Lambda = (tau_plus - 1/2) (tau_minus - 1/2)
            // Choose such that HWBB walls are always in the right place.
            // TODO: make this a configurable parameter.
            const distribn_t Lambda = 3.0 / 16.0;

            const distribn_t tau_plus = lbmParams->GetTau();
            const distribn_t omega_plus = lbmParams->GetOmega();
            const distribn_t tau_minus = 0.5 + Lambda / (tau_plus - 0.5);
            const distribn_t omega_minus =  -1.0 / tau_minus;

            if constexpr (HasZero) {
                // Special case the null velocity.
                hydroVars.SetFPostCollision(iZero,
                                            hydroVars.f[iZero] + omega_plus * hydroVars.f_neq.f[iZero]);
            }

            // Now deal with the non-zero
            for (auto [i, iBar]: directionPairs)
            {
                distribn_t sym = 0.5 * omega_plus * (hydroVars.f_neq.f[i] + hydroVars.f_neq.f[iBar]);
                distribn_t asym = 0.5 * omega_minus * (hydroVars.f_neq.f[i] - hydroVars.f_neq.f[iBar]);
                hydroVars.SetFPostCollision(i, hydroVars.f[i] + sym + asym);
                hydroVars.SetFPostCollision(iBar, hydroVars.f[iBar] + sym - asym);
            }
        }
    };
}

#endif /* HEMELB_LB_KERNELS_TRT_H */
