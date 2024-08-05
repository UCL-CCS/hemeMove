// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_IO_READERS_XDRREADER_H
#define HEMELB_IO_READERS_XDRREADER_H

#include <cstdint>
#include <span>
#include <string>

#include "Exception.h"
#include "io/XdrSerialisation.h"
#include "util/Vector3D.h"

namespace hemelb::io
{
    // Base class to read XDR data. This does the XDR related
    // stuff. Derived classes must implement GetPosition and get_bytes
    // member functions to allow it to work.
    class XdrReader
    {
    public:
        // Virtual destructor.
        virtual ~XdrReader() = default;

        // Main function for reading the next bit of the stream.
        //
        // Stores result of deserialising in the argument supplied and
        // indicates success with return value.
        template<class T>
        bool read(T& val) {
            constexpr auto n = xdr::xdr_serialised_size<T>();
            auto buf = get_bytes(n);
            xdr::xdr_deserialise(val, buf);
            return true;
        }

        // Like xdr_opaque
        template<std::size_t N>
        bool read(std::span<std::byte, N>& bytes);

        template <typename T>
        bool read(util::Vector3D<T>& vec) {
            for (int i = 0; i < 3; ++i)
                read(vec[i]);
            return true;
        }

        // Return result of deserialising - will throw on error.
        template <class T>
        T read() {
            T ans;
            if (!read<T>(ans)) {
                throw Exception() << "Error reading type from XDR";
            }
            return ans;
        }

        // Get the position in the stream.
        virtual unsigned GetPosition() = 0;

    protected:
        // Get some bytes from the underlying storage
        virtual const std::byte* get_bytes(size_t n) = 0;
    };

    // Specialisation for strings
    template<>
    inline bool XdrReader::read(std::string& val) {
        uint32_t len = 0;
        read(len);
        // p == padded
        uint32_t plen = 4 * ((len - 1)/4 + 1);
        auto pstr = reinterpret_cast<const char*>(get_bytes(plen));
        val.assign(pstr, len);
        return true;
    }
    // Like xdr opaque
    template<std::size_t N>
    inline bool XdrReader::read(std::span<std::byte, N>& bytes) {
        const auto len = bytes.size();
        const auto req_nwords = (len - 1)/4 + 1;
        const auto space = req_nwords * 4;
        auto data = get_bytes(space);
        std::copy(data, data + len, bytes.begin());
        return true;
    }
}

#endif // HEMELB_IO_READERS_XDRREADER_H
