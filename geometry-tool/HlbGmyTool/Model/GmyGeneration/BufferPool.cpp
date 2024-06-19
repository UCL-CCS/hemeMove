// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "BufferPool.h"
#include <cstddef>

namespace hemelb::gmytool::gmy {

BufferPool::BufferPool(unsigned int bSize) : size(bSize) {}

BufferPool::~BufferPool() {
  // Clear the stack
  while (!this->unused.empty()) {
    delete[] this->unused.top();
    this->unused.pop();
  }
}

char* BufferPool::New() {
  // If the stack is empty, create a new array, otherwise pop an array
  if (this->unused.empty()) {
    return new char[this->size];
  } else {
    char* ans = this->unused.top();
    this->unused.pop();
    return ans;
  }
}

void BufferPool::Free(char* buf) {
  // If the buffer is NULL, skip
  if (buf == NULL)
    return;
  // If we have fewer than 10, add this one to the unused, otherwise delete it
  if (this->unused.size() < 10) {
    this->unused.push(buf);
  } else {
    delete[] buf;
  }
}

// Return the size of buffers handled.
unsigned int BufferPool::GetSize() const {
  return this->size;
}

}  // namespace hemelb::gmytool::gmy
