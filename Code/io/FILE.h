// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_IO_FILE_H
#define HEMELB_IO_FILE_H

#include <cstdio>
#include <filesystem>
#include <memory>

namespace hemelb::io {
    // A thin RAII-style owning wrapper around a C standard library FILE object.
    // Note this is NOT MPI-aware, so mainly for use in tests.
    class FILE {
        struct Closer {
            void operator()(std::FILE *);
        };

        using ClosingPtr = std::unique_ptr<std::FILE, Closer>;
        ClosingPtr _handle;
    public:
        static FILE open(std::filesystem::path const &p, std::string_view mode);

        // Direct file read - calls std::fread on the underlying file handle
        std::size_t read(void *buffer, std::size_t size, std::size_t count);

        // Direct file read - calls std::fwrite on the underlying file handle
        std::size_t write(void const *buffer, std::size_t size, std::size_t count);

        long tell() const;

        int seek(long offset, long origin);

        void clearerr();
    };
}
#endif