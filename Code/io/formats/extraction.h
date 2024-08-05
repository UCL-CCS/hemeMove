// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_IO_FORMATS_EXTRACTION_H
#define HEMELB_IO_FORMATS_EXTRACTION_H

namespace hemelb::io::formats
{
    struct extraction {
        // Magic number to identify extraction data files.
        // ASCII for 'xtr' + EOF
        static constexpr std::uint32_t MagicNumber = 0x78747204;

        // The version number of the file format.
        static constexpr std::uint32_t VersionNumber = 6;

        // The length of the main header. Made up of:
        // uint - HemeLbMagicNumber
        // uint - ExtractionMagicNumber
        // uint - Format version number
        // (Parameters for conversion to physical units)
        // double - dx (metres)
        // double - dt (seconds)
        // double - dm (kg)
        // double x 3 - Origin x,y,z components (metres)
        // double - reference pressure (Pa)
        // uhyper - Total number of sites
        // uint - Field count
        // uint - Length of the field header that follows
        static constexpr std::size_t MainHeaderLength = 84;

        // Each field header is:
        // string - name
        // uint32 - number of elements
        // uint32 - type code
        // uint32 - number of offsets (valid values are {0, 1, n elem})
        // type[n offsets] - offsets (n offsets items of type implied above)
        // type - scale lattice to physical units (0 => no scaling)
        enum class TypeCode : std::uint32_t {
            FLOAT,
            DOUBLE,
            INT32,
            UINT32,
            INT64,
            UINT64,
        };

        // Compute the length of data written by XDR for a given string.
        static inline std::size_t GetStoredLengthOfString(std::string const &str) {
            // XDR pads up to the nearest multiple of four bytes and also
            // stores the length of the string.
            auto len_bytes = str.length();
            auto len_words = (len_bytes - 1U) / 4U + 1U;
            return (1 + len_words) * 4U;
        }

        // Compute the length (in bytes) of a field header
        static inline std::size_t GetFieldHeaderLength(std::string const &name, std::uint32_t noff, TypeCode tc) {
            std::size_t len = GetStoredLengthOfString(name);
            len += 4U;  // number of elements
            len += 4U;  // type code
            len += 4U;  // number of offsets
            std::size_t elemsize = (tc == TypeCode::FLOAT || tc == TypeCode::INT32 || tc == TypeCode::UINT32) ? 4U : 8U;
            len += elemsize * noff; // offsets
            len += elemsize; // scale
            return len;
        }
    };
}
#endif // HEMELB_IO_FORMATS_EXTRACTION_H
