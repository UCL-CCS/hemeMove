// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELBSETUPTOOL_BLOCKWRITER_H
#define HEMELBSETUPTOOL_BLOCKWRITER_H

#include <stddef.h>
#include <string>
#include "io/writers/XdrMemWriter.h"

class GeometryWriter;
class BufferPool;
/*
 * Extension of a hemelb::io::XdrWriter that notes how many fluid sites, in how
 * much space, have been written. It then pushes this to the GeometryWriter's
 * headerEncoder and writes to the GeometryWriter's body
 */

class BlockWriter {
 public:
  BlockWriter(BufferPool* bp);
  void Reset();

  ~BlockWriter();

  void IncrementFluidSitesCount();

  void Finish();
  void Write(GeometryWriter& gw);

  // Overload << to delegate to the XdrMemWriter
  template <typename T>
  BlockWriter& operator<<(T const& value) {
    (*this->writer) << value;
    return *this;
  }

 protected:
  std::byte* buffer;
  hemelb::io::XdrMemWriter* writer;
  BufferPool* bufferPool;
  unsigned nFluidSites;
  unsigned CompressedBlockLength;
  unsigned UncompressedBlockLength;
  bool IsFinished;
};

#endif  // HEMELBSETUPTOOL_BLOCKWRITER_H
