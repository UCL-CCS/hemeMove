// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "io/hexfloat.h"
#include "io/ensure_hexfloat.h"
#include <sstream>
#include "Exception.h"
#include "log/Logger.h"

namespace hemelb::io {

    namespace {
        constexpr char _tolower(char c) {
            if (c >= 'A' && c <= 'Z')
                return c + 32;
            return c;
        }
        constexpr bool _is_dec_digit(char c) {
            return c >= '0' && c <= '9';
        }
        constexpr bool _is_hex_digit(char c) {
            return _is_dec_digit(c) || (c >= 'a' && c <= 'f');
        }
    }

    // We have a horrible state machine to determine which
    // characters are allowed, as have decimal and hexadecimal
    // floats to handle.
    //
    // Consider letters as lowercase.
    //
    // Dec = decimal, Hex = hexadecimal
    // Man == mantissa, exp = exponent
    // Pre and post refer to the decimal point
    //
    // Any condition not handled is end of number
    //
    // Start: +/- => ManSign, 0 => UnknownMan, 1-9 => DecManPre
    //
    // ManSign: 0 => UnknownMan, 1-9 => DecManPre
    //
    // UnknownMan: x => HexManPre, 0-9 => DecManPre, . => DecDP
    //
    // DecManPre: 0-9 => DecManPre, . => DecDP, e => DecExpStart
    //
    // DecDP: 0-9 => DecManPost, e => DecExpStart
    // DecManPost: 0-9 => DecManPost, e => DecExpStart
    // NOTE: since these have the same transitions, we can merge
    //
    // DecExpStart: +- => DecExp, 0-9 => DecExp
    //
    // DecExp: 0-9 => DecExp
    //
    // HexManPre: 0-f => HexManPre, . => HexDP, p => DecExpStart
    //
    // HexDP: 0-f => HexManPost, p => DecExpStart
    //
    // HexManPost: 0-f => HexManPost, p => DecExpStart
    //
    // NOTE: yes, hex floats have decimal exponents!

    bool FloatCharAccumulator::add_char(char c) {
        char lc = _tolower(c);
        using enum State;

        switch (state) {

            case Start:
                if (c == '+' || c == '-') {
                    state = ManSign;
                } else if (c == '0') {
                    state = UnknownMan;
                } else if (c > '0' && c <= '9') {
                    state = DecManPre;
                } else {
                    state = Error;
                }
                break;

            case ManSign:
                if (c == '0') {
                    state = UnknownMan;
                } else if (c > '0' && c <= '9') {
                    state = DecManPre;
                } else {
                    state = Error;
                }
                break;

            case UnknownMan:
                if (lc == 'x') {
                    state = HexManPre;
                } else if (_is_dec_digit(c)) {
                    state = DecManPre;
                } else if (c == '.') {
                    state = DecManPost; // Would be DecDP but same as this.
                } else {
                    state = Error;
                }
                break;

            case DecManPre:
                if (_is_dec_digit(c)) {
                    state = DecManPre;
                } else if (c == '.') {
                    state = DecManPost; // Would be DecDP but same as this.
                } else if (lc == 'e') {
                    state = DecExpStart;
                } else {
                    state = Error;
                }
                break;

            case DecManPost:
                if (_is_dec_digit(c)) {
                    state = DecManPost;
                } else if (lc == 'e') {
                    state = DecExpStart;
                } else {
                    state = Error;
                }
                break;

            case DecExpStart:
                if (c == '+' || c == '-' || _is_dec_digit(c)) {
                    state = DecExp;
                } else {
                    state = Error;
                }
                break;

            case DecExp:
                if (_is_dec_digit(c)) {
                    state = DecExp;
                } else {
                    state = Error;
                }
                break;

            case HexManPre:
                if (_is_hex_digit(lc)) {
                    state = HexManPre;
                } else if (c == '.') {
                    state = HexManPost; // Would be HexDP but same as this.
                } else if (lc == 'p') {
                    state = DecExpStart;
                } else {
                    state = Error;
                }
                break;

            case HexManPost:
                if (_is_hex_digit(lc)) {
                    state = HexManPost;
                } else if (lc == 'p') {
                    state = DecExpStart;
                } else {
                    state = Error;
                }
                break;

            default:
                throw (Exception() << "Unknown case in hexfloat character accumulation");
        }

        if (state == Error) {
            return true;
        } else {
            buf.push_back(c);
            return false;
        }
    }

    bool GlobalHexFloatLocale::CurrentCanParseHexFloats() {
        // 0.3 can't be exactly represented as FP
        double expected = 0.3;
        auto rep = "0x1.3333333333333p-2";
        std::istringstream s(rep);
        double read = 0;
        s >> read;
        return read == expected;
    }

    struct GlobalHexFloatLocale::Impl {
        std::locale original;
    };
    GlobalHexFloatLocale::GlobalHexFloatLocale() : impl(std::make_unique<Impl>()) {
        if (!CurrentCanParseHexFloats()) {
            log::Logger::Log<log::Info, log::Singleton>(
                    "Default locale does not parse hexadecimal floats, installing custom global facets");
            auto hexloc = std::locale(
                    std::locale(impl->original, new reliable_hexfloat_numeric_facet<char>),
                    new reliable_hexfloat_numeric_facet<wchar_t>
            );
            std::locale::global(hexloc);
        }
    }

    GlobalHexFloatLocale::~GlobalHexFloatLocale() {
        std::locale::global(impl->original);
    }
}