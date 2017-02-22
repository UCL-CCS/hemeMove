
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_EXTRACTION_PROPERTYACTOR_H
#define HEMELB_EXTRACTION_PROPERTYACTOR_H

#include "extraction/PropertyWriter.h"
#include "io/PathManager.h"
#include "lb/MacroscopicPropertyCache.h"
#include "lb/SimulationState.h"
#include "timestep/Actor.h"

namespace hemelb
{
  namespace extraction
  {
    class PropertyActor : public timestep::Actor
    {
      public:
        /**
         * Constructor, gets the class ready for reading.
         * @param simulationState
         * @param propertyOutputs
         * @param dataSource
         * @return
         */
        PropertyActor(const lb::SimulationState& simulationState,
                      const std::vector<PropertyOutputFile*>& propertyOutputs,
                      IterableDataSource& dataSource,
                      reporting::Timers& timers,
                      comm::Communicator::ConstPtr ioComms);

        ~PropertyActor();

        /**
         * Set which properties will be required this iteration.
         * @param propertyCache
         */
        void SetRequiredProperties(lb::MacroscopicPropertyCache& propertyCache);

      inline virtual void BeginAll() {
      }
      inline virtual void Begin()  {
      }
      inline virtual void Receive()  {
      }
      inline virtual void PreSend()  {
      }
      inline virtual void Send()  {
      }
      inline virtual void PreWait()  {
      }
      inline virtual void Wait() {
      }
      inline virtual void End()  {
      }
      
      // Override the iterated actor end of iteration method to perform writing.
      virtual void EndAll();
      
      private:
        const lb::SimulationState& simulationState;
        PropertyWriter* propertyWriter;
        reporting::Timers& timers;
    };
  }
}

#endif /* HEMELB_EXTRACTION_PROPERTYACTOR_H */
