// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_IO_WRITERS_XDRMEMWRITER_H
#define HEMELB_IO_WRITERS_XDRMEMWRITER_H

#include "io/writers/XdrWriter.h"

namespace hemelb::io
{
    // Simple wrapper to allow construction from pointer-length
    // arguments.
    class XdrMemWriter : public XdrMetaWriter<std::byte*> {
    public:
        XdrMemWriter(std::byte* dataBuffer, unsigned int dataLength);
    };
}
#endif // HEMELB_IO_WRITERS_XDRMEMWRITER_H
