// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "io/FILE.h"
#include <cstring>
#include "Exception.h"

namespace hemelb::io {
    void FILE::Closer::operator()(std::FILE *fh) {
        std::fclose(fh);
    }

    FILE FILE::open(std::filesystem::path const &p, std::string_view mode) {
        FILE ans;
        ans._handle = ClosingPtr(std::fopen(p.c_str(), mode.data()));
        if (!ans._handle.get())
            throw Exception() << "Error opening file '" << p << "' with reason: " << std::strerror(errno);
        return ans;
    }

    // Direct file read - calls std::fread on the underlying file handle
    std::size_t FILE::read(void *buffer, std::size_t size, std::size_t count) {
        return std::fread(buffer, size, count, _handle.get());
    }

    // Direct file read - calls std::fwrite on the underlying file handle
    std::size_t FILE::write(void const *buffer, std::size_t size, std::size_t count) {
        return std::fwrite(buffer, size, count, _handle.get());
    }

    long FILE::tell() const {
        return std::ftell(_handle.get());
    }

    int FILE::seek(long offset, long origin) {
        return std::fseek(_handle.get(), offset, origin);
    }

    void FILE::clearerr() {
        std::clearerr(_handle.get());
    }
}