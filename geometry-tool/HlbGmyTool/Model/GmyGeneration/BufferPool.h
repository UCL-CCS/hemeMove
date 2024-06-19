// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HLBGMYTOOL_GMY_BUFFERPOOL_H
#define HLBGMYTOOL_GMY_BUFFERPOOL_H
#include <stack>

namespace hemelb::gmytool::gmy {

// Allocates and frees or reuses buffers of a given size.
class BufferPool {
 public:
  // C'tor- argument is the size of buffers to work with.
  BufferPool(unsigned int);
  ~BufferPool();
  // Returns an uninitialized buffer
  char* New();
  // Return a buffer to the pool or free it
  void Free(char*);
  // Returns the size of buffers managed by the pool
  unsigned int GetSize() const;

 private:
  unsigned int size;
  std::stack<char*> unused;
};

}  // namespace hemelb::gmytool::gmy
#endif
