// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_IO_FORMATS_GEOMETRY_H
#define HEMELB_IO_FORMATS_GEOMETRY_H

#include <array>

#include "io/formats/formats.h"
#include "util/Vector3D.h"

namespace hemelb::io::formats
{
    // Class that contains information necessary to interpret
    // a HemeLB geometry file (*.gmy).
    //
    // The file is described in Doc/dev/formats/GmyReadResult.md
    //
    // Types and values should match those used in the file.
    struct geometry
    {
        // Magic number to identify geometry files.
        // ASCII for 'gmy', then EOF
        // Combined magic number is:
        // hex    68 6c 62 21 67 6d 79 04
        // ascii:  h  l  b  !  g  m  y EOF
        static constexpr std::uint32_t MagicNumber = 0x676d7904;

        // Version number for geometry format.
        static constexpr std::uint32_t VersionNumber = 4;

        // Type codes permitted for sites
        enum class SiteType : std::uint32_t {
            SOLID = 0,
            FLUID = 1
        };

        // Type codes for the sort of boundaries intersected by a link.
        enum class CutType : std::uint32_t
        {
            NONE = 0,
            WALL = 1,
            INLET = 2,
            OUTLET = 3
        };

        // Type codes defining wall normal availability
        enum class WallNormalAvailability : std::uint32_t
        {
            NOT_AVAILABLE = 0,
            AVAILABLE = 1
        };

        // Number of displacements in the neighbourhood
        static constexpr size_t NumberOfDisplacements = 26;

        // The length of the preamble for the geometry file:
        //  * 1 uint for the HemeLB magic number
        //  * 1 uint for the geometry magic number
        //  * 1 uint for the version
        //  * 3 uints for the problem dimensions in blocks
        //  * 1 uint for the number of sites along one block side
        //  * 1 uint, value 0 to pad to 32 bytes
        //
        //  * 8 uints = 8 * 4  = 32
        static constexpr size_t PreambleLength = 32;

        // The length of a single header record (i.e. you have one of
        // these per block):
        //  * 1 uint for number of fluid sites
        //  * 1 uint for compressed data size in bytes
        //  * 1 uint for uncompressed data size in bytes
        static constexpr size_t HeaderRecordLength = 12;

        // The maximum possible length of a single site's data.
        //  * 1 uint for the type
        //  * then NumberOfDisplacements:
        //    * 1 uint for the cut type
        //    * 1 uint for the inlet/outlet ID
        //    * 1 float for the cut distance
        //  * 1 uint for the wall normal availability
        //  * 3 floats for the wall normal
        static constexpr size_t MaxFluidSiteRecordLength = 4 + geometry::NumberOfDisplacements * (4 + 4 + 4) + 4 + 3 * 4;

        // Maximum length of a solid site's data. (In fact, this is
        // THE length of a solid site's data.)
        static constexpr size_t MaxSolidSiteRecordLength = 4;

        // Compute the maximum possible length of a single block's data.
        // @param blockSideLength
        // @return maximum block record length in bytes
        static constexpr size_t GetMaxBlockRecordLength(size_t blockSideLength)
        {
            return blockSideLength * blockSideLength * blockSideLength * geometry::MaxFluidSiteRecordLength;
        }

        // Compute the maximum length of a block's data, given the
        // number of fluid sites contained within it.
        // @param blockSideLength
        // @param nFluidSites
        // @return max length in bytes
        static constexpr size_t GetMaxBlockRecordLength(size_t blockSideLength, size_t nFluidSites)
        {
            size_t nSolidSites = blockSideLength * blockSideLength * blockSideLength - nFluidSites;
            return (nFluidSites * geometry::MaxFluidSiteRecordLength + nSolidSites * geometry::MaxSolidSiteRecordLength);
        }

        // Aliases for displacement vectors to neighbouring sites.
        using Displacement = util::Vector3D<int>;
        using DisplacementArray = std::array<Displacement, NumberOfDisplacements>;

        // The ordering of vectors making up the 3D Moore
        // neighbourhood of a site that is used in the file format.
        static constexpr DisplacementArray Neighbourhood{{
	    Displacement(-1, -1, -1),
	    Displacement(-1, -1, 0),
	    Displacement(-1, -1, +1),

	    Displacement(-1, 0, -1),
	    Displacement(-1, 0, 0),
	    Displacement(-1, 0, +1),

	    Displacement(-1, +1, -1),
	    Displacement(-1, +1, 0),
	    Displacement(-1, +1, +1),

	    Displacement(0, -1, -1),
	    Displacement(0, -1, 0),
	    Displacement(0, -1, +1),

	    Displacement(0, 0, -1),
	    // Displacement( 0, 0, 0) - missing as there is never anything to say about the link to self
	    Displacement(0, 0, +1),

	    Displacement(0, +1, -1),
	    Displacement(0, +1, 0),
	    Displacement(0, +1, +1),

	    Displacement(+1, -1, -1),
	    Displacement(+1, -1, 0),
	    Displacement(+1, -1, +1),

	    Displacement(+1, 0, -1),
	    Displacement(+1, 0, 0),
	    Displacement(+1, 0, +1),

	    Displacement(+1, +1, -1),
	    Displacement(+1, +1, 0),
	    Displacement(+1, +1, +1)
	  }};
    };
}
#endif // HEMELB_IO_FORMATS_GEOMETRY_H
