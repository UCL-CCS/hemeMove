// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "Index.h"

namespace hemelb {

namespace gmytool {

namespace {
void DefaultHandlerFunction(int direction) {
  throw IndexError();
}
}  // namespace

}  // namespace gmytool

namespace util {
Vector3DBase::HandlerFunction* Vector3DBase::handler =
    gmytool::DefaultHandlerFunction;
}
}  // namespace hemelb
