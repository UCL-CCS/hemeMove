#ifndef HEMELB_EXTRACTION_PROPERTYWRITER_H
#define HEMELB_EXTRACTION_PROPERTYWRITER_H

#include "extraction/LocalPropertyOutput.h"
#include "extraction/PropertyOutputFile.h"
#include "mpiInclude.h"

namespace hemelb
{
  namespace extraction
  {
    class PropertyWriter
    {
      public:
        /**
         * Constructor, takes a vector of the output files to create.
         * @param propertyOutputs
         * @return
         */
        PropertyWriter(IterableDataSource& dataSource, const std::vector<PropertyOutputFile>& propertyOutputs);

        /**
         * Destructor; deallocates memory used to store property info.
         * @return
         */
        ~PropertyWriter();

        /**
         * Writes each of the property output files, if appropriate for the passed iteration number.
         *
         * An iterationNumber of 0 will write all files.
         * @param iterationNumber
         */
        void Write(unsigned long iterationNumber) const;

      private:
        /**
         * Holds sufficient information to output property information from this core.
         */
        std::vector<LocalPropertyOutput*> localPropertyOutputs;
    };
  }
}

#endif /* HEMELB_EXTRACTION_PROPERTYWRITER_H */
