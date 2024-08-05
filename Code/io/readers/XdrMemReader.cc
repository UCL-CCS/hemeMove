// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "io/readers/XdrMemReader.h"
#include "hassert.h"

namespace hemelb::io
{
  // Constructor to create an Xdr object based on a memory buffer
  XdrMemReader::XdrMemReader(const std::byte* buf, unsigned int dataLength)
    : start(buf), current(buf), len(dataLength)
  {
  }

  XdrMemReader::XdrMemReader(const std::vector<std::byte>& dataVec)
    : start(dataVec.data()), current(start), len(dataVec.size())
  {
  }

  unsigned XdrMemReader::GetPosition() {
    return current - start;
  }

  const std::byte* XdrMemReader::get_bytes(size_t n) {
    HASSERT(GetPosition() + n <= len);
    auto ans = current;
    current += n;
    return ans;
  }

}
