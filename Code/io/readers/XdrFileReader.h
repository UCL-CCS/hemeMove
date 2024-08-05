// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_IO_READERS_XDRFILEREADER_H
#define HEMELB_IO_READERS_XDRFILEREADER_H

#include <cstdio>
#include <filesystem>
#include <memory>
#include <vector>
#include "io/readers/XdrReader.h"
#include "io/FILE.h"

namespace hemelb::io
{
  // Deserialise from a file given by path.
  class XdrFileReader : public XdrReader {
  public:
    XdrFileReader(const std::filesystem::path& fn);
    ~XdrFileReader() override = default;
    unsigned GetPosition() override;
  protected:
    const std::byte* get_bytes(size_t n) override;
  private:
    FILE fh;
    // Buffer for holding read data
    std::vector<std::byte> buf;
  };

}

#endif  // HEMELB_IO_READERS_XDRFILEREADER_H
