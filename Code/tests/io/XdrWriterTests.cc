// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <type_traits>

#include <catch2/catch.hpp>

#include "io/writers/XdrWriter.h"
// This header is generated by xdr_gen.py
#include "tests/io/xdr_test_data.h"

namespace hemelb::tests
{

    template<typename T>
    void TestBasic(const std::vector<T>& values, const std::vector<std::byte>& expected_buffer) {
      // Figure out sizes and counts
      constexpr auto sz = sizeof(T);
      static_assert(sz > 0 && sz <= 8,
		    "Only works on 8--64 bit types");
      constexpr auto type_bits = sz*8;
      // XDR works on 32 bit words
      constexpr auto coded_words = (sz - 1)/4 + 1;
      constexpr auto coded_bytes = coded_words * 4;
      // Need a value for each bit being on + zero
      constexpr auto n_vals = type_bits + 1;
      // Buffers for the encoded data
      constexpr auto buf_size = n_vals * coded_bytes;
      std::byte our_buf[buf_size];
      // Fill with binary ones to trigger failure if we're not
      // writing enough zeros
      std::fill(our_buf, our_buf+buf_size, std::byte{std::numeric_limits<unsigned char>::max()});

      auto our_coder = io::MakeXdrWriter(our_buf, our_buf + buf_size);
	  
      REQUIRE(n_vals == values.size());
      for (auto& val: values) {
	our_coder << val;
      }
	  
      REQUIRE(buf_size == expected_buffer.size());
      for (auto i = 0U; i < buf_size; ++i) {
	REQUIRE(expected_buffer[i] == our_buf[i]);
      }
    }


    TEMPLATE_TEST_CASE("XdrWriter works for integers", "", int32_t, uint32_t, int64_t, uint64_t) {
      
      // Flip every bit in the range in turn and check it encodes to
      // the same as libc's xdr_encode. Works only for 1- or 2-word
      // encoded size types.
      static_assert(std::is_integral<TestType>::value, "only works on (u)ints");
      auto&& values = test_data<TestType>::unpacked();
      auto&& expected_buffer = test_data<TestType>::packed();
      TestBasic(values, expected_buffer);
    }

    TEMPLATE_TEST_CASE("XdrWriter works for floating types", "", float, double) {
      static_assert(std::is_floating_point<TestType>::value, "Floats only please!");

      auto&& values = test_data<TestType>::unpacked();
      auto&& expected_buffer = test_data<TestType>::packed();
	  
      constexpr bool is_double = sizeof(TestType) == 8;
      REQUIRE(TestType{0.0} == values[0]);
      constexpr auto sign_bit = sizeof(TestType)*8 - 1;
      REQUIRE(-TestType{0.0} == values[sign_bit + 1]);

      REQUIRE(std::numeric_limits<TestType>::denorm_min() == values[1]);

      constexpr size_t mantissa_bits = is_double ? 52 : 23;
      REQUIRE(std::numeric_limits<TestType>::min() == values[mantissa_bits+1]);
      TestBasic(values, expected_buffer);
    }

    TEST_CASE("XdrWriter works for strings") {
        using UPB = std::unique_ptr<std::byte[]>;
        auto make_ones = [](size_t n) -> UPB {
            // Fill with binary ones to trigger failure if we're not
            // writing enough zeros
            auto ans = UPB(new std::byte[n]);
            std::fill(ans.get(), ans.get() + n, std::byte{std::numeric_limits<unsigned char>::max()});
            return ans;
        };

        // XDR strings are serialised as length, data (0-padded to a word)
        auto coded_length = [](const std::string& s) {
            return s.empty() ? 4U : 4U * (2U + (s.size() - 1U) / 4U);
        };

        auto&& values = test_data<std::string>::unpacked();
        auto&& expected_buffer = test_data<std::string>::packed();

        std::size_t total_length = 0;
        std::for_each(
                values.begin(), values.end(),
                [&](const std::string& v) {
                    total_length += coded_length(v);
                });
        REQUIRE(expected_buffer.size() == total_length);

        auto our_buf = make_ones(total_length);
        auto our_coder = io::MakeXdrWriter(our_buf.get(), our_buf.get() + total_length);

        for (auto& v: values) {
            our_coder << v;
        }
        for (auto i = 0U; i < total_length; ++i) {
            CAPTURE(i);
            REQUIRE(expected_buffer[i] == our_buf[i]);
        }
    }
}
