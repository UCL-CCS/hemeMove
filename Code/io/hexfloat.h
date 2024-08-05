// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_IO_HEXFLOAT_H
#define HEMELB_IO_HEXFLOAT_H

#include <algorithm>
#include <concepts>
#include <locale>
#include <string>

// Parsing of hexfloats is not implemented in libstdc++, see
// https://gcc.gnu.org/bugzilla/show_bug.cgi?id=81122
// Hence, all this...
namespace hemelb::io {

    struct FloatCharAccumulator {
        enum class State {
            Start,
            ManSign,
            UnknownMan,
            DecManPre,
            // DecDP, // See note in implementation.
            DecManPost,
            DecExpStart,
            DecExp,
            HexManPre,
            //HexDP, // See note in implementation.
            HexManPost,
            Error
        };

        std::string buf;
        State state = State::Start;
        // Accumulate the character, return if done
        bool add_char(char c);
    };

    template <
            class CharT,
            class InputIt = std::istreambuf_iterator<CharT>
    >
    class reliable_hexfloat_numeric_facet : public std::num_get<CharT, InputIt> {
    public:
        using base_type = std::num_get<CharT, InputIt>;
        using char_type = CharT;

        using base_type::base_type;

    protected:
        ~reliable_hexfloat_numeric_facet() override = default;

        InputIt do_get( InputIt in, InputIt end, std::ios_base& str,
                        std::ios_base::iostate& err, float& v ) const override {
            return convert_float("%g", std::strtof, in, end, str, err, v);
        }
        InputIt do_get( InputIt in, InputIt end, std::ios_base& str,
                        std::ios_base::iostate& err, double& v ) const override {
            return convert_float("%lg", std::strtod, in, end, str, err, v);
        }
        InputIt do_get( InputIt in, InputIt end, std::ios_base& str,
                        std::ios_base::iostate& err, long double& v ) const override {
            return convert_float("%Lg", std::strtold, in, end, str, err, v);
        }

        template <std::floating_point T>
        InputIt convert_float(
                char const* specifier, T (conv_func)(const char*, char**),
                InputIt& in, InputIt const& end,
                std::ios_base& str, std::ios_base::iostate& err,
                T& v
        ) const {
            // Stage 1 done by setting parameters to this function

            // Stage 2
            auto loc = str.getloc();

            // atoms for a float
            static const char src[] = "0123456789abcdefpxABCDEFPX+-";
            constexpr auto N = sizeof(src) - 1;
            char_type atoms[sizeof(src)];
            auto const& char_facet = std::use_facet<std::ctype<char_type>>(str.getloc());
            char_facet.widen(src, src + sizeof(src), atoms);
            auto const& np_facet = std::use_facet<std::numpunct<char_type>>(loc);
            const auto dp = np_facet.decimal_point();
            const auto ks = np_facet.thousands_sep();
            const auto gl = np_facet.grouping().length();

            bool have_dp = false;

            FloatCharAccumulator acca;
            for (; in != end; ++in) {
                char_type ct = *in;

                char c = src[std::find(atoms, atoms + N, ct) - atoms];
                if (ct == dp) {
                    c = '.';
                    have_dp = true;
                }

                if (ct == ks && gl != 0) {
                    if (!have_dp) {
                        // Should remember the position and discard the char (but not checking thousands)
                        continue;
                    } else {
                        // End stage 2
                        break;
                    }
                }

                if (c == 0)
                    // Not allowed (find got end == '\0'), end S2
                    break;

                if (acca.add_char(c))
                    break;
            }

            // Stage 3 - accumulated chars in acca.buf
            char* conv_end;
            v = conv_func(acca.buf.data(), &conv_end);
            if (std::distance(acca.buf.data(), conv_end) != std::ssize(acca.buf)) {
                v = 0;
                err = std::ios_base::failbit;
            }
            // Should check digit grouping, but no
            if (in == end)
                err |= std::ios_base::eofbit;
            return in;
        }
    };

    template <typename CharT, typename Traits>
    void HexImbue(std::basic_istream<CharT, Traits>& stream) {
        auto hexloc = std::locale(
                stream.getloc(),
                new reliable_hexfloat_numeric_facet<CharT>
        );
        stream.imbue(hexloc);
    }
}

#endif