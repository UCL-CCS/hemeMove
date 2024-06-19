// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HLBGMYTOOL_GMY_GEOMETRYWRITER_H
#define HLBGMYTOOL_GMY_GEOMETRYWRITER_H

#include <cstdio>
#include <string>

#include "Index.h"

#include "io/writers/xdr/XdrWriter.h"

namespace hemelb::gmytool::gmy {

using io::writers::xdr::XdrWriter;

class BlockWriter;
class BufferPool;

class GeometryWriter {
 public:
  GeometryWriter(const std::string& OutputGeometryFile,
                 int BlockSize,
                 Index BlockCounts);

  ~GeometryWriter();

  void Close();
  BlockWriter* StartNextBlock();

 protected:
  std::string OutputGeometryFile;
  int BlockSize;
  Index BlockCounts;

  int headerStart;
  XdrWriter* headerEncoder;
  unsigned int headerBufferLength;
  char* headerBuffer;

  int bodyStart;
  FILE* bodyFile;
  BufferPool* BlockBufferPool;
  friend class BlockWriter;
};

}  // namespace hemelb::gmytool::gmy
#endif  // HLBGMYTOOL_GMY_GEOMETRYWRITER_H
