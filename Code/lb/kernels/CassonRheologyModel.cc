// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "lb/kernels/CassonRheologyModel.h"
#include "util/numerical.h"
#include <cmath>

namespace hemelb::lb
{
    double CassonRheologyModel::CalculateViscosityForShearRate(const double &iShearRate,
                                                               const distribn_t &iDensity) const
    {
        double k0_k1_gamma = K0 + K1 * std::sqrt(iShearRate);
        double eta = (k0_k1_gamma * k0_k1_gamma) / iShearRate;

        // In the Casson rheology model, viscosity tends to infinity as shear rate goes to zero. Bound it.
        return std::min(eta, CASSON_MAX_VISCOSITY);
    }
}
