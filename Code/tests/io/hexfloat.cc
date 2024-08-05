// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <sstream>
#include <catch2/catch.hpp>

#include "io/hexfloat.h"
#include "util/Vector3D.h"

namespace hemelb::tests
{
    using io::reliable_hexfloat_numeric_facet;

    namespace {
        template <typename... Ts>
        auto MakeStream(Ts&&... args) {
            auto stream = std::istringstream(std::forward<Ts>(args)...);
            io::HexImbue(stream);
            return stream;
        }
        template <typename T>
        auto Parse(std::string const& s) {
            auto N = s.length();
            auto istream = MakeStream(s);
            T ans;
            istream >> ans;
            auto pos = istream.tellg();
            // EOF => -1
            if (istream.eof())
                pos = N;
            return std::make_pair(ans, pos);
        }
    }
    template <typename T>
    void CheckValid(T expected, char const* repr) {
        std::string s(repr);
        auto [val, N] = Parse<T>(s);
        REQUIRE(N == std::ssize(s));
        REQUIRE(val == expected);
    }

// Use the macros for stringification
#define CheckValidDouble(x) CheckValid<double>(x, #x)
// Not space so can have the text arg include the parens
#define CheckValidVec(x) \
    CheckValid<util::Vector3D<double>>(util::Vector3D x, #x)

    TEST_CASE("hexfloat") {
        SECTION("Successful parses to double") {
            CheckValidDouble(0.0);
            CheckValidDouble(1.0);
            CheckValidDouble(+0.3);
            CheckValidDouble(-4.5214e+4);
            CheckValidDouble(+157.47521E-14);
            CheckValidDouble(0x0.001p0);
            CheckValidDouble(-0xa.234fP+12);
            CheckValidDouble(+0xDEADBEEFp-1);
        }

        SECTION("Failing parses to double") {
            auto [val, N] = Parse<double>("0.1a");
            REQUIRE(val == 0.1);
            REQUIRE(N == 3);
        }

        SECTION("Vectors of doubles") {
            CheckValidVec((+0.0,0.0e0,0x0p0));
            CheckValidVec((1.0,-4.129e+3,-0Xa.145b69P+2));
        }
    }
}