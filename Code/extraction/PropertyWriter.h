// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_EXTRACTION_PROPERTYWRITER_H
#define HEMELB_EXTRACTION_PROPERTYWRITER_H

#include "extraction/LocalPropertyOutput.h"
#include "extraction/PropertyOutputFile.h"
#include "net/mpi.h"

namespace hemelb::extraction
{
    class PropertyWriter
    {
      public:
        /**
         * Constructor, takes a vector of the output files to create.
         * @param propertyOutputs
         * @return
         */
        PropertyWriter(std::shared_ptr<IterableDataSource> dataSource,
                       const std::vector<PropertyOutputFile>& propertyOutputs,
                       net::IOCommunicator const& ioComms);

        /**
         * Writes each of the property output files, if appropriate for the passed iteration number.
         *
         * An iterationNumber of 0 will write all files.
         * @param iterationNumber
         */
        void Write(unsigned long iterationNumber, unsigned long totalSteps);

        /**
         * Returns a vector of all the LocalPropertyOutputs.
         * @return
         */
        const std::vector<LocalPropertyOutput>& GetPropertyOutputs() const;

      private:
        /**
         * Holds sufficient information to output property information from this core.
         */
        std::vector<LocalPropertyOutput> localPropertyOutputs;
    };
}

#endif /* HEMELB_EXTRACTION_PROPERTYWRITER_H */
