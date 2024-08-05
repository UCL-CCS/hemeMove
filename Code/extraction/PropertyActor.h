// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_EXTRACTION_PROPERTYACTOR_H
#define HEMELB_EXTRACTION_PROPERTYACTOR_H

#include "extraction/PropertyWriter.h"
#include "lb/MacroscopicPropertyCache.h"
#include "lb/SimulationState.h"
#include "net/IteratedAction.h"
#include "reporting/timers_fwd.h"

namespace hemelb::extraction
{
    class PropertyActor : public net::IteratedAction
    {
      public:
        /**
         * Constructor, gets the class ready for reading.
         * @param simulationState
         * @param propertyOutputs
         * @param dataSource
         * @return
         */
        PropertyActor(std::shared_ptr<lb::SimulationState const> simulationState,
                      const std::vector<PropertyOutputFile>& propertyOutputs,
                      std::shared_ptr<IterableDataSource> dataSource, reporting::Timers& timers,
                      const net::IOCommunicator& ioComms);

        ~PropertyActor() override;

        /**
         * Set which properties will be required this iteration.
         * @param propertyCache
         */
        void SetRequiredProperties(lb::MacroscopicPropertyCache& propertyCache);

        /**
         * Override the iterated actor end of iteration method to perform writing.
         */
        void EndIteration() override;

      private:
        std::shared_ptr<lb::SimulationState const> simulationState;
        std::unique_ptr<PropertyWriter> propertyWriter;
        reporting::Timers& timers;
    };
}

#endif /* HEMELB_EXTRACTION_PROPERTYACTOR_H */
