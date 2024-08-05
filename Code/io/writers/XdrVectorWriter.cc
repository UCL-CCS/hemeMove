// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "io/writers/XdrVectorWriter.h"

namespace hemelb::io
{
	// Use the magic protected constructor.
	// First function obj makes a std::vector<char>
	// Second one creates a std::back_insert_iterator from it
	XdrVectorWriter::XdrVectorWriter() : base([](){
						    return byte_vec{};
						  },
						  [](byte_vec& v){
						    return std::back_inserter(v);
						  })
	{
	}
	// Delegate to above then do the reserve
	XdrVectorWriter::XdrVectorWriter(size_t n) : XdrVectorWriter()
	{
	  res.reserve(n);
	}

	const byte_vec& XdrVectorWriter::GetBuf() const
	{
	  return res;
	}

}
