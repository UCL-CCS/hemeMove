# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

if(NOT HAVE_XDRUINTXX_T)
  #Without the standard names for the xdr functions, create aliases for the existing ones.
  add_definitions(-Dxdr_uint16_t=xdr_u_int16_t)
  add_definitions(-Dxdr_uint32_t=xdr_u_int32_t)
  add_definitions(-Dxdr_uint64_t=xdr_u_int64_t)
endif()

add_library(hemelb_io OBJECT
  FILE.cc
  writers/AsciiFileWriter.cc writers/AsciiStreamWriter.cc
  readers/XdrFileReader.cc
  readers/XdrMemReader.cc
  writers/XdrWriter.cc
  writers/XdrMemWriter.cc
  writers/XdrFileWriter.cc
  writers/XdrVectorWriter.cc
  writers/NullWriter.cc
  writers/Writer.cc
  xml.cc
  Checkpointer.cc
  TimePattern.cc
  hexfloat.cc
)

target_link_libraries(hemelb_io PUBLIC
  hemelb_util
  )
target_link_libraries(hemelb_io PRIVATE
  tinyxml2::tinyxml2
  )
