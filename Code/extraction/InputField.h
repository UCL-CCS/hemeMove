// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_EXTRACTION_INPUTFIELD_H
#define HEMELB_EXTRACTION_INPUTFIELD_H

#include <cstdint>
#include <string>
#include "io/formats/extraction.h"

namespace hemelb
{
  namespace extraction
  {
    struct InputField
    {
      std::string name;
      std::uint32_t numberOfElements;
      std::uint32_t typecode;
      std::uint32_t numberOfOffsets;
      // The type of offsets varies depending on typecode, but for now
      // we are only dealing with double precision data (the
      // distributions) without an offset, so just omit.
      //
      // std::vector<T> offsets;

      // Also have a scale factor from written lattice units -> physical.
      // Again, we are dealing with distributions in lattice units so omit.
      // T scale;
    };
  }
}

#endif /* HEMELB_EXTRACTION_INPUTFIELD_H */
