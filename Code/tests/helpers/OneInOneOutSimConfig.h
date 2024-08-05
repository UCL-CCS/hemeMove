// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_TESTS_HELPERS_ONEINONEOUTSIMCONFIG_H
#define HEMELB_TESTS_HELPERS_ONEINONEOUTSIMCONFIG_H

#include <string>

#include "io/xml.h"
#include "configuration/SimConfigReader.h"
#include "configuration/SimBuilder.h"

namespace hemelb::tests::helpers
{
    // TODO: Figure out what this is supposed to be.
    class OneInOneOutSimConfigReader : public configuration::SimConfigReader
    {
    public:
        OneInOneOutSimConfigReader();
        void CheckIoletMatchesCMake(const io::xml::Element& ioletEl,
                                    std::string_view requiredBC) const override;

        [[nodiscard]] configuration::SimConfig Read() const override;
        [[nodiscard]] configuration::GlobalSimInfo
        DoIOForSimulation(const io::xml::Element simEl) const override;
        [[nodiscard]] configuration::ICConfig DoIOForInitialConditions(io::xml::Element parent) const override;
        [[nodiscard]] std::vector<configuration::IoletConfig>
        DoIOForInOutlets(configuration::GlobalSimInfo const& sim_info, const io::xml::Element xmlNode) const override;

    };

}

#endif
