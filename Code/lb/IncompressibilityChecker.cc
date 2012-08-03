// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#include "lb/lattices/D3Q15.h"
#include "lb/IncompressibilityChecker.hpp"

namespace hemelb
{
  namespace lb
  {
    // Explicit instantiation
    template class IncompressibilityChecker<net::PhasedBroadcastRegular<> > ;
  }
}
