// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "io/writers/XdrFileWriter.h"

namespace hemelb::io
{

    BufferedFileWriter BufferedFileWriter::open(const std::filesystem::path &p) {
        BufferedFileWriter ans;
        ans.f = FILE::open(p, "w");
        ans.buf.resize(1024);
        ans.pos = 0;
        return ans;
    }

    void BufferedFileWriter::put(const std::byte &c) {
        if (pos + 1 == buf.size())
            flush();
        buf[pos] = c;
        ++pos;
    }

    void BufferedFileWriter::flush() {
        if (pos) {
            f.write(buf.data(), 1, pos);
            pos = 0;
        }
    }

    long BufferedFileWriter::tell() const {
        return f.tell() + pos;
    }

    BufferedFileWriter::~BufferedFileWriter() noexcept {
        flush();
    }

    BufferedFileWriterIter::BufferedFileWriterIter(BufferedFileWriter* w) : file{w} {}

    BufferedFileWriterIter& BufferedFileWriterIter::operator*() {
        return *this;
    }
    BufferedFileWriterIter& BufferedFileWriterIter::operator++() {
        return *this;
    }
    BufferedFileWriterIter& BufferedFileWriterIter::operator++(int) {
        return *this;
    }
    BufferedFileWriterIter& BufferedFileWriterIter::operator=(std::byte const& b) {
        file->put(b);
        return *this;
    }
    // Use base class special constructor
    // Functor 1 opens the file
    // Functor 2 makes an iterator that writes to it
    XdrFileWriter::XdrFileWriter(const std::filesystem::path& fileName) :
            base([&fileName]() {
                     return BufferedFileWriter::open(fileName);
                 },
                 [](BufferedFileWriter& fb){
                     return BufferedFileWriterIter(&fb);
                 })
    {
    }

}
