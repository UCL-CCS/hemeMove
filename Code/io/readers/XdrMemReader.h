// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_IO_READERS_XDRMEMREADER_H
#define HEMELB_IO_READERS_XDRMEMREADER_H

#include <vector>

#include "io/readers/XdrReader.h"

namespace hemelb::io
{

  // Deserialise from some memory.
  //
  // NB: this class does not OWN the memory you pass to its
  // constructor. Caller is responsible for ensuring its lifetime
  // exceeds that of the object.
  class XdrMemReader : public XdrReader
  {
  public:
    XdrMemReader(const std::byte* dataBuffer, unsigned int dataLength);
    XdrMemReader(const std::vector<std::byte>& dataVec);
    ~XdrMemReader() override = default;
    unsigned GetPosition() override;
  protected:
    const std::byte* get_bytes(size_t n) override;
  private:
    const std::byte* start;
    const std::byte* current;
    size_t len;
  };
}

#endif  // HEMELB_IO_READERS_XDRMEMREADER_H
