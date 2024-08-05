// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <cstdio>

#include <zlib.h>

#include "BlockWriter.h"
#include "BufferPool.h"
#include "GenerationError.h"
#include "GeometryWriter.h"
#include "Neighbours.h"

BlockWriter::BlockWriter(BufferPool* bp)
    : writer(NULL), buffer(NULL), bufferPool(bp) {
  this->Reset();
}

void BlockWriter::Reset() {
  this->nFluidSites = 0;
  this->CompressedBlockLength = 0;
  this->UncompressedBlockLength = 0;

  this->IsFinished = false;

  this->bufferPool->Free(this->buffer);
  this->buffer = this->bufferPool->New();

  delete this->writer;
  this->writer =
      new hemelb::io::XdrMemWriter(this->buffer, this->bufferPool->GetSize());
}

BlockWriter::~BlockWriter() {
  delete this->writer;
  this->bufferPool->Free(this->buffer);
}

void BlockWriter::IncrementFluidSitesCount() {
  this->nFluidSites++;
}

void BlockWriter::Finish() {
  this->CompressedBlockLength = 0;
  this->UncompressedBlockLength = 0;

  if (this->nFluidSites > 0) {
    int ret;  // zlib return code

    // How much data to compress?
    this->UncompressedBlockLength = this->writer->getCurrentStreamPosition();

    // Set up our compressor
    z_stream stream;
    stream.zalloc = Z_NULL;
    stream.zfree = Z_NULL;
    stream.opaque = Z_NULL;
    // Max compression
    ret = deflateInit(&stream, 9);
    if (ret != Z_OK)
      throw GenerationErrorMessage("Cannot init zlib structures");

    //		// Compute upper bound for how much space we'll need.
    //		int maxDeflatedLength = deflateBound(&stream,
    //				this->UncompressedBlockLength);
    auto compressedBuffer = this->bufferPool->New();

    // Set input. The XDR buffer has to be char but zlib only works with
    // unsigned char. Just cast for now...
    stream.next_in = reinterpret_cast<unsigned char*>(this->buffer);
    stream.avail_in = this->UncompressedBlockLength;
    // Set output
    stream.next_out = reinterpret_cast<unsigned char*>(compressedBuffer);
    stream.avail_out = this->bufferPool->GetSize();

    // Deflate. This should be it, if not their was an error.
    ret = deflate(&stream, Z_FINISH);
    if (ret != Z_STREAM_END)
      throw GenerationErrorMessage("Error compressing buffer");

    // How much space did we actually use?
    this->CompressedBlockLength =
        reinterpret_cast<std::byte*>(stream.next_out) - compressedBuffer;

    // Tell zlib to clean up.
    ret = deflateEnd(&stream);
    if (ret != Z_OK)
      throw GenerationErrorMessage("Cannot free zlib structures");

    std::swap(this->buffer, compressedBuffer);
    this->bufferPool->Free(compressedBuffer);
  } else {
    this->bufferPool->Free(this->buffer);
    this->buffer = NULL;
  }
  delete this->writer;
  this->writer = NULL;

  this->IsFinished = true;
}

void BlockWriter::Write(GeometryWriter& gw) {
  if (this->nFluidSites > 0) {
    if (this->buffer == NULL)
      throw GenerationErrorMessage("Cannot write NULL buffer");

    std::fwrite(this->buffer, 1, this->CompressedBlockLength, gw.bodyFile);
  }
  *(gw.headerEncoder) << this->nFluidSites << this->CompressedBlockLength
                      << this->UncompressedBlockLength;
}
