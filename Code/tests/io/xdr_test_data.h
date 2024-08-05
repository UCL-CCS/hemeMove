// -*- mode: c++; -*-
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_TESTS_IO_XDR_TEST_DATA_H
#define HEMELB_TESTS_IO_XDR_TEST_DATA_H
#include <vector>

namespace hemelb::tests
{
    template <typename T>
    struct test_data {
      static const std::vector<T>& unpacked();
      static const std::vector<std::byte>& packed();
    };
}
#endif
