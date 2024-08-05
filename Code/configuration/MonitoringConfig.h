// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_CONFIGURATION_MONITORINGCONFIG_H
#define HEMELB_CONFIGURATION_MONITORINGCONFIG_H

#include "extraction/PropertyOutputFile.h"

namespace hemelb::configuration
{

    // Bundles together various configuration parameters concerning simulation monitoring
    struct MonitoringConfig
    {
        bool doConvergenceCheck = false; ///< Whether to turn on the convergence check or not
        extraction::source::Type convergenceVariable; ///< Macroscopic variable used to check for convergence
        double convergenceReferenceValue; ///< Reference value used to normalise an absolute error (making it relative)
        double convergenceRelativeTolerance = 0.0; ///< Convergence check relative tolerance
        bool convergenceTerminate = false; ///< Whether to terminate a converged run or not
        bool doIncompressibilityCheck = false; ///< Whether to turn on the IncompressibilityChecker or not
    };
}

#endif
