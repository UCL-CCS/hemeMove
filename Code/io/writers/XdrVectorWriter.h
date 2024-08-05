// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_IO_WRITERS_XDRVECTORWRITER_H
#define HEMELB_IO_WRITERS_XDRVECTORWRITER_H

#include <vector>

#include "io/writers/XdrWriter.h"

namespace hemelb::io
{
    using byte_vec = std::vector<std::byte>;
    using byte_vec_appender = std::back_insert_iterator<byte_vec>;

    // XDR encoder that will put its buffer in a vector so you
    // don't have to worry about managing memory. The encoder owns
    // the memory but you can get a reference to the data with
    // GetBuf()
    class XdrVectorWriter : public XdrMetaWriter<byte_vec_appender, byte_vec>
    {
        using base = XdrMetaWriter<byte_vec_appender, byte_vec>;
    public:
        // Start with no memory reserved
        XdrVectorWriter();
        // Reserve n bytes of memory for storage
        XdrVectorWriter(size_t n);

        const byte_vec& GetBuf() const;
    };

}

#endif
