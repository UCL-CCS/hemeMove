// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#include <cassert>
#include <cstring>

#include "io/readers/XdrFileReader.h"

namespace hemelb::io
{
  // Constructor to create an Xdr object based on a file.
  XdrFileReader::XdrFileReader(const std::filesystem::path& fn) :
    fh{FILE::open(fn, "r")} {
  }

  unsigned XdrFileReader::GetPosition() {
    return fh.tell();
  }

  const std::byte* XdrFileReader::get_bytes(size_t n) {
    buf.clear();
    buf.resize(n);
    auto nread = fh.read(buf.data(), 1, n);
    if (nread != n)
        throw Exception() << "Could not read " << n << " bytes, instead got " << nread;
    return buf.data();
  }

}
