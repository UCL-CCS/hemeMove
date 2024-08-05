// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_IO_FORMATS_COLLOIDS_H
#define HEMELB_IO_FORMATS_COLLOIDS_H

#include "io/formats/formats.h"

namespace hemelb
{
  namespace io
  {
    namespace formats
    {
      namespace colloids
      {
        /* The format comprises a header and then a body.
         * Header is made up of (hex file position, type, description)
         * 00   uint       HemeLB magic number (see formats.h)
         * 04   uint       Colloids magic number (see below)
         * 08   uint       Version number
         * 12   uint       Offset of the body start into the file (bytes)
         * 16   uint       Flag for simulation stability
         * 20   double     Voxel size in physical units (m)
         * 28   3 x dbl    Origin of voxel indexing in physical units (m)
         * 52   3 x uint   Bounding box minimum indices
         * 64   3 x uint   Bounding box maximum indices
         * 76   uint       Number of fluid voxels
         * Header length = 80 bytes
         *
         * Body consists of one record per voxel
         * 3 x uint     Voxel index
         * 1 x float   Pressure (Pa)
         * 3 x float   Velocity (m/s)
         * 1 x float   von Mises stress (Pa). Note this is NaN for non-wall sites
         * Record length = 32 bytes
         *
         */

        enum
        {
          /* Identify colloid files
           * ASCII for 'col', then EOF
           * Combined magic number is
           * hex    68 6c 62 21 63 6f 6c 04
           * ascii:  h  l  b  !  c  o  l EOF
           */
          MagicNumber = 0x636f6c04
        };
        enum
        {
          VersionNumber = 1
        };
        enum
        {
          MagicLength = 12
        };
        enum
        {
          HeaderLength = 24
        };

        enum
        {
          RecordLength = 56
        };
      }
    }
  }

}
#endif // HEMELB_IO_FORMATS_COLLOIDS_H
