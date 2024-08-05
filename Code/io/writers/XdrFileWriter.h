// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_IO_WRITERS_XDRFILEWRITER_H
#define HEMELB_IO_WRITERS_XDRFILEWRITER_H

#include <string>
#include <io/FILE.h>
#include "io/writers/XdrWriter.h"

namespace hemelb::io
{
    // Simple buffering of writing unformatted bytes to a file
    class BufferedFileWriter {
        io::FILE f;
        std::vector<std::byte> buf;
        std::size_t pos = 0;

        BufferedFileWriter() = default;

    public:
        static BufferedFileWriter open(std::filesystem::path const &p);

        // No copy
        BufferedFileWriter(BufferedFileWriter const&) = delete;
        BufferedFileWriter& operator=(BufferedFileWriter const&) = delete;
        // Move is OK
        BufferedFileWriter(BufferedFileWriter&&) = default;
        BufferedFileWriter& operator=(BufferedFileWriter&&) = default;
        // Destruction flushes and closes
        ~BufferedFileWriter() noexcept;

        void put(std::byte const& c);
        void flush();

        long tell() const;
    };

    // Basic output iterator for use by XdrMetaWriter
    class BufferedFileWriterIter {
        BufferedFileWriter* file;
    public:
        using difference_type = std::ptrdiff_t;

        BufferedFileWriterIter(BufferedFileWriter*);

        BufferedFileWriterIter& operator*();
        BufferedFileWriterIter& operator++();
        BufferedFileWriterIter& operator++(int);
        BufferedFileWriterIter& operator=(std::byte const& b);
    };

    // Wrapper that owns a filehandle and uses that as the sink for data
    class XdrFileWriter : public XdrMetaWriter<BufferedFileWriterIter, BufferedFileWriter>
    {
        using base = XdrMetaWriter<BufferedFileWriterIter, BufferedFileWriter>;
    public:
        explicit XdrFileWriter(const std::filesystem::path& fileName);
    };

}
#endif // HEMELB_IO_WRITERS_XDRFILEWRITER_H
