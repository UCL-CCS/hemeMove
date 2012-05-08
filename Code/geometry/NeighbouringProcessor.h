#ifndef HEMELB_GEOMETRY_NEIGHBOURINGPROCESSOR_H
#define HEMELB_GEOMETRY_NEIGHBOURINGPROCESSOR_H

#include "units.h"

namespace hemelb
{
  namespace geometry
  {
    /**
     * Encapsulates data about a processor which has at least 1 site which
     * neighbours at least one site local to this core.
     */
    struct NeighbouringProcessor
    {
      public:
        // Rank of the neighbouring processor.
        proc_t Rank;

        // The number of distributions shared between this neighbour and the current processor.
        site_t SharedFCount;

        // Index on this processor of the first distribution shared between this
        // neighbour and the current processor.
        site_t FirstSharedF;
    };
  }
}

#endif /* HEMELB_GEOMETRY_NEIGHBOURINGPROCESSOR_H */
