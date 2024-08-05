// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_TESTS_HELPERS_SIMCONFBUILDHELP_H
#define HEMELB_TESTS_HELPERS_SIMCONFBUILDHELP_H

#include <memory>
#include "configuration/SimConfig.h"
#include "configuration/SimBuilder.h"

namespace hemelb::tests {

    class UninitSimBuilder : public configuration::SimBuilder {
    public:
        inline UninitSimBuilder(configuration::SimConfig const &conf, std::shared_ptr<util::UnitConverter> conv)
                : configuration::SimBuilder{conf, false} {
            unit_converter = std::move(conv);
        }
    };
}
#endif // ONCE
