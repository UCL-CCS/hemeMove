// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "tests/io/xdr_test_data.h"

#include <cstdint>
#include <cstring>
#include <string>

namespace hemelb::tests
{

    namespace {
        float binflt(std::uint32_t x) {
            float ans;
            std::memcpy(&ans, &x, 4);
            return ans;
        }
        double bindbl(std::uint64_t x) {
            double ans;
            std::memcpy(&ans, &x, 8);
            return ans;
        }
    }
    template<>
    const std::vector<std::int32_t>& test_data<std::int32_t>::unpacked() {
      static const std::vector<std::int32_t> unpacked = {
        0x0,
        0x1,
        0x2,
        0x4,
        0x8,
        0x10,
        0x20,
        0x40,
        0x80,
        0x100,
        0x200,
        0x400,
        0x800,
        0x1000,
        0x2000,
        0x4000,
        0x8000,
        0x10000,
        0x20000,
        0x40000,
        0x80000,
        0x100000,
        0x200000,
        0x400000,
        0x800000,
        0x1000000,
        0x2000000,
        0x4000000,
        0x8000000,
        0x10000000,
        0x20000000,
        0x40000000,
        static_cast<std::int32_t>(0x80000000)
      };
      return unpacked;
    }
    
    template<>
    const std::vector<std::byte>& test_data<std::int32_t>::packed() {
      static const std::vector<std::byte> packed = {std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{1}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{2}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{4}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{8}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{16}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{32}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{64}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{128}, std::byte{0}, std::byte{0}, std::byte{1}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{2}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{4}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{8}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{16}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{32}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{64}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{128}, std::byte{0}, std::byte{0}, std::byte{1}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{2}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{4}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{8}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{16}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{32}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{64}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{128}, std::byte{0}, std::byte{0}, std::byte{1}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{2}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{4}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{8}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{16}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{32}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{64}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{128}, std::byte{0}, std::byte{0}, std::byte{0}};
      return packed;
    };
template<>
    const std::vector<std::int64_t>& test_data<std::int64_t>::unpacked() {
      static const std::vector<std::int64_t> unpacked = {
        0x0,
        0x1,
        0x2,
        0x4,
        0x8,
        0x10,
        0x20,
        0x40,
        0x80,
        0x100,
        0x200,
        0x400,
        0x800,
        0x1000,
        0x2000,
        0x4000,
        0x8000,
        0x10000,
        0x20000,
        0x40000,
        0x80000,
        0x100000,
        0x200000,
        0x400000,
        0x800000,
        0x1000000,
        0x2000000,
        0x4000000,
        0x8000000,
        0x10000000,
        0x20000000,
        0x40000000,
        0x80000000,
        0x100000000,
        0x200000000,
        0x400000000,
        0x800000000,
        0x1000000000,
        0x2000000000,
        0x4000000000,
        0x8000000000,
        0x10000000000,
        0x20000000000,
        0x40000000000,
        0x80000000000,
        0x100000000000,
        0x200000000000,
        0x400000000000,
        0x800000000000,
        0x1000000000000,
        0x2000000000000,
        0x4000000000000,
        0x8000000000000,
        0x10000000000000,
        0x20000000000000,
        0x40000000000000,
        0x80000000000000,
        0x100000000000000,
        0x200000000000000,
        0x400000000000000,
        0x800000000000000,
        0x1000000000000000,
        0x2000000000000000,
        0x4000000000000000,
        static_cast<std::int64_t>(0x8000000000000000)
      };
      return unpacked;
    }
    
    template<>
    const std::vector<std::byte>& test_data<std::int64_t>::packed() {
      static const std::vector<std::byte> packed = {std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{1}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{2}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{4}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{8}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{16}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{32}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{64}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{128}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{1}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{2}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{4}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{8}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{16}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{32}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{64}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{128}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{1}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{2}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{4}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{8}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{16}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{32}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{64}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{128}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{1}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{2}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{4}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{8}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{16}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{32}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{64}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{128}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{1}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{2}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{4}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{8}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{16}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{32}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{64}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{128}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{1}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{2}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{4}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{8}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{16}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{32}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{64}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{128}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{1}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{2}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{4}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{8}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{16}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{32}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{64}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{128}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{1}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{2}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{4}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{8}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{16}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{32}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{64}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{128}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}};
      return packed;
    };
template<>
    const std::vector<std::uint32_t>& test_data<std::uint32_t>::unpacked() {
      static const std::vector<std::uint32_t> unpacked = {
        0x0,
        0x1,
        0x2,
        0x4,
        0x8,
        0x10,
        0x20,
        0x40,
        0x80,
        0x100,
        0x200,
        0x400,
        0x800,
        0x1000,
        0x2000,
        0x4000,
        0x8000,
        0x10000,
        0x20000,
        0x40000,
        0x80000,
        0x100000,
        0x200000,
        0x400000,
        0x800000,
        0x1000000,
        0x2000000,
        0x4000000,
        0x8000000,
        0x10000000,
        0x20000000,
        0x40000000,
        0x80000000
      };
      return unpacked;
    }
    
    template<>
    const std::vector<std::byte>& test_data<std::uint32_t>::packed() {
      static const std::vector<std::byte> packed = {std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{1}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{2}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{4}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{8}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{16}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{32}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{64}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{128}, std::byte{0}, std::byte{0}, std::byte{1}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{2}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{4}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{8}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{16}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{32}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{64}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{128}, std::byte{0}, std::byte{0}, std::byte{1}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{2}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{4}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{8}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{16}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{32}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{64}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{128}, std::byte{0}, std::byte{0}, std::byte{1}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{2}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{4}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{8}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{16}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{32}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{64}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{128}, std::byte{0}, std::byte{0}, std::byte{0}};
      return packed;
    };
template<>
    const std::vector<std::uint64_t>& test_data<std::uint64_t>::unpacked() {
      static const std::vector<std::uint64_t> unpacked = {
        0x0,
        0x1,
        0x2,
        0x4,
        0x8,
        0x10,
        0x20,
        0x40,
        0x80,
        0x100,
        0x200,
        0x400,
        0x800,
        0x1000,
        0x2000,
        0x4000,
        0x8000,
        0x10000,
        0x20000,
        0x40000,
        0x80000,
        0x100000,
        0x200000,
        0x400000,
        0x800000,
        0x1000000,
        0x2000000,
        0x4000000,
        0x8000000,
        0x10000000,
        0x20000000,
        0x40000000,
        0x80000000,
        0x100000000,
        0x200000000,
        0x400000000,
        0x800000000,
        0x1000000000,
        0x2000000000,
        0x4000000000,
        0x8000000000,
        0x10000000000,
        0x20000000000,
        0x40000000000,
        0x80000000000,
        0x100000000000,
        0x200000000000,
        0x400000000000,
        0x800000000000,
        0x1000000000000,
        0x2000000000000,
        0x4000000000000,
        0x8000000000000,
        0x10000000000000,
        0x20000000000000,
        0x40000000000000,
        0x80000000000000,
        0x100000000000000,
        0x200000000000000,
        0x400000000000000,
        0x800000000000000,
        0x1000000000000000,
        0x2000000000000000,
        0x4000000000000000,
        0x8000000000000000
      };
      return unpacked;
    }
    
    template<>
    const std::vector<std::byte>& test_data<std::uint64_t>::packed() {
      static const std::vector<std::byte> packed = {std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{1}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{2}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{4}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{8}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{16}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{32}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{64}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{128}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{1}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{2}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{4}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{8}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{16}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{32}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{64}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{128}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{1}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{2}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{4}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{8}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{16}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{32}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{64}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{128}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{1}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{2}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{4}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{8}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{16}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{32}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{64}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{128}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{1}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{2}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{4}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{8}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{16}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{32}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{64}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{128}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{1}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{2}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{4}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{8}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{16}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{32}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{64}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{128}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{1}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{2}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{4}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{8}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{16}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{32}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{64}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{128}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{1}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{2}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{4}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{8}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{16}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{32}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{64}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{128}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}};
      return packed;
    };
template<>
    const std::vector<float>& test_data<float>::unpacked() {
      static const std::vector<float> unpacked = {
        binflt(0x0),
        binflt(0x1),
        binflt(0x2),
        binflt(0x4),
        binflt(0x8),
        binflt(0x10),
        binflt(0x20),
        binflt(0x40),
        binflt(0x80),
        binflt(0x100),
        binflt(0x200),
        binflt(0x400),
        binflt(0x800),
        binflt(0x1000),
        binflt(0x2000),
        binflt(0x4000),
        binflt(0x8000),
        binflt(0x10000),
        binflt(0x20000),
        binflt(0x40000),
        binflt(0x80000),
        binflt(0x100000),
        binflt(0x200000),
        binflt(0x400000),
        binflt(0x800000),
        binflt(0x1000000),
        binflt(0x2000000),
        binflt(0x4000000),
        binflt(0x8000000),
        binflt(0x10000000),
        binflt(0x20000000),
        binflt(0x40000000),
        binflt(0x80000000)
      };
      return unpacked;
    }
    
    template<>
    const std::vector<std::byte>& test_data<float>::packed() {
      static const std::vector<std::byte> packed = {std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{1}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{2}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{4}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{8}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{16}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{32}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{64}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{128}, std::byte{0}, std::byte{0}, std::byte{1}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{2}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{4}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{8}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{16}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{32}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{64}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{128}, std::byte{0}, std::byte{0}, std::byte{1}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{2}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{4}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{8}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{16}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{32}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{64}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{128}, std::byte{0}, std::byte{0}, std::byte{1}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{2}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{4}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{8}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{16}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{32}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{64}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{128}, std::byte{0}, std::byte{0}, std::byte{0}};
      return packed;
    };
template<>
    const std::vector<double>& test_data<double>::unpacked() {
      static const std::vector<double> unpacked = {
        bindbl(0x0),
        bindbl(0x1),
        bindbl(0x2),
        bindbl(0x4),
        bindbl(0x8),
        bindbl(0x10),
        bindbl(0x20),
        bindbl(0x40),
        bindbl(0x80),
        bindbl(0x100),
        bindbl(0x200),
        bindbl(0x400),
        bindbl(0x800),
        bindbl(0x1000),
        bindbl(0x2000),
        bindbl(0x4000),
        bindbl(0x8000),
        bindbl(0x10000),
        bindbl(0x20000),
        bindbl(0x40000),
        bindbl(0x80000),
        bindbl(0x100000),
        bindbl(0x200000),
        bindbl(0x400000),
        bindbl(0x800000),
        bindbl(0x1000000),
        bindbl(0x2000000),
        bindbl(0x4000000),
        bindbl(0x8000000),
        bindbl(0x10000000),
        bindbl(0x20000000),
        bindbl(0x40000000),
        bindbl(0x80000000),
        bindbl(0x100000000),
        bindbl(0x200000000),
        bindbl(0x400000000),
        bindbl(0x800000000),
        bindbl(0x1000000000),
        bindbl(0x2000000000),
        bindbl(0x4000000000),
        bindbl(0x8000000000),
        bindbl(0x10000000000),
        bindbl(0x20000000000),
        bindbl(0x40000000000),
        bindbl(0x80000000000),
        bindbl(0x100000000000),
        bindbl(0x200000000000),
        bindbl(0x400000000000),
        bindbl(0x800000000000),
        bindbl(0x1000000000000),
        bindbl(0x2000000000000),
        bindbl(0x4000000000000),
        bindbl(0x8000000000000),
        bindbl(0x10000000000000),
        bindbl(0x20000000000000),
        bindbl(0x40000000000000),
        bindbl(0x80000000000000),
        bindbl(0x100000000000000),
        bindbl(0x200000000000000),
        bindbl(0x400000000000000),
        bindbl(0x800000000000000),
        bindbl(0x1000000000000000),
        bindbl(0x2000000000000000),
        bindbl(0x4000000000000000),
        bindbl(0x8000000000000000)
      };
      return unpacked;
    }
    
    template<>
    const std::vector<std::byte>& test_data<double>::packed() {
      static const std::vector<std::byte> packed = {std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{1}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{2}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{4}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{8}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{16}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{32}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{64}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{128}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{1}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{2}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{4}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{8}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{16}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{32}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{64}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{128}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{1}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{2}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{4}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{8}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{16}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{32}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{64}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{128}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{1}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{2}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{4}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{8}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{16}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{32}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{64}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{128}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{1}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{2}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{4}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{8}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{16}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{32}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{64}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{128}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{1}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{2}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{4}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{8}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{16}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{32}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{64}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{128}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{1}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{2}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{4}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{8}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{16}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{32}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{64}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{128}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{1}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{2}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{4}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{8}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{16}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{32}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{64}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{128}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}};
      return packed;
    };
template<>
    const std::vector<std::string>& test_data<std::string>::unpacked() {
      static const std::vector<std::string> unpacked = {
        "",
        "a",
        "test",
        "This is a slightly longer string whose length is a multiple of four.",
        "This is a slightly longer string whose length isn't a multiple of four."
      };
      return unpacked;
    }
    
    template<>
    const std::vector<std::byte>& test_data<std::string>::packed() {
      static const std::vector<std::byte> packed = {std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{1}, std::byte{97}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{4}, std::byte{116}, std::byte{101}, std::byte{115}, std::byte{116}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{68}, std::byte{84}, std::byte{104}, std::byte{105}, std::byte{115}, std::byte{32}, std::byte{105}, std::byte{115}, std::byte{32}, std::byte{97}, std::byte{32}, std::byte{115}, std::byte{108}, std::byte{105}, std::byte{103}, std::byte{104}, std::byte{116}, std::byte{108}, std::byte{121}, std::byte{32}, std::byte{108}, std::byte{111}, std::byte{110}, std::byte{103}, std::byte{101}, std::byte{114}, std::byte{32}, std::byte{115}, std::byte{116}, std::byte{114}, std::byte{105}, std::byte{110}, std::byte{103}, std::byte{32}, std::byte{119}, std::byte{104}, std::byte{111}, std::byte{115}, std::byte{101}, std::byte{32}, std::byte{108}, std::byte{101}, std::byte{110}, std::byte{103}, std::byte{116}, std::byte{104}, std::byte{32}, std::byte{105}, std::byte{115}, std::byte{32}, std::byte{97}, std::byte{32}, std::byte{109}, std::byte{117}, std::byte{108}, std::byte{116}, std::byte{105}, std::byte{112}, std::byte{108}, std::byte{101}, std::byte{32}, std::byte{111}, std::byte{102}, std::byte{32}, std::byte{102}, std::byte{111}, std::byte{117}, std::byte{114}, std::byte{46}, std::byte{0}, std::byte{0}, std::byte{0}, std::byte{71}, std::byte{84}, std::byte{104}, std::byte{105}, std::byte{115}, std::byte{32}, std::byte{105}, std::byte{115}, std::byte{32}, std::byte{97}, std::byte{32}, std::byte{115}, std::byte{108}, std::byte{105}, std::byte{103}, std::byte{104}, std::byte{116}, std::byte{108}, std::byte{121}, std::byte{32}, std::byte{108}, std::byte{111}, std::byte{110}, std::byte{103}, std::byte{101}, std::byte{114}, std::byte{32}, std::byte{115}, std::byte{116}, std::byte{114}, std::byte{105}, std::byte{110}, std::byte{103}, std::byte{32}, std::byte{119}, std::byte{104}, std::byte{111}, std::byte{115}, std::byte{101}, std::byte{32}, std::byte{108}, std::byte{101}, std::byte{110}, std::byte{103}, std::byte{116}, std::byte{104}, std::byte{32}, std::byte{105}, std::byte{115}, std::byte{110}, std::byte{39}, std::byte{116}, std::byte{32}, std::byte{97}, std::byte{32}, std::byte{109}, std::byte{117}, std::byte{108}, std::byte{116}, std::byte{105}, std::byte{112}, std::byte{108}, std::byte{101}, std::byte{32}, std::byte{111}, std::byte{102}, std::byte{32}, std::byte{102}, std::byte{111}, std::byte{117}, std::byte{114}, std::byte{46}, std::byte{0}};
      return packed;
    };

}
