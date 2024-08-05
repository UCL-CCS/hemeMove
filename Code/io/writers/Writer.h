// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_IO_WRITERS_WRITER_H
#define HEMELB_IO_WRITERS_WRITER_H

#include <cstdint>
#include <span>
#include <string>

#include "util/Vector3D.h"

namespace hemelb::io
{
      class Writer
      {
        public:

          enum Separator
          {
            eol
          };

          // Special version for eol, using function overloading
          Writer& operator<<(enum Separator const & value);

          // Overload << to write basic types and any necessary separators
          template<typename T>
          Writer& operator<<(T const & value)
          {
            _write(value);
            writeFieldSeparator();
            return *this;
          }

          template <typename T>
          Writer& operator<<(util::Vector3D<T> const& vec)
          {
              return *this << vec.x() << vec.y() << vec.z();
          }

          // Function to get the current position of writing in the stream.
          virtual unsigned int getCurrentStreamPosition() const = 0;

          virtual ~Writer() = default;
        protected:
          Writer() = default;

          // Functions for formatting control
          virtual void writeFieldSeparator() = 0;
          virtual void writeRecordSeparator() = 0;

          // Methods to simply write (no separators) which are virtual and
          // hence must be overriden.
          virtual void _write(std::int16_t const& intToWrite) = 0;
          virtual void _write(std::uint16_t const& uIntToWrite) = 0;
          virtual void _write(std::int32_t const& intToWrite) = 0;
          virtual void _write(std::uint32_t const& uIntToWrite) = 0;
          virtual void _write(std::int64_t const& intToWrite) = 0;
          virtual void _write(std::uint64_t const& uIntToWrite) = 0;

          virtual void _write(double const& doubleToWrite) = 0;
          virtual void _write(float const& floatToWrite) = 0;

          virtual void _write(const std::string& floatToWrite) = 0;
          virtual void _write(std::span<const std::byte>) = 0;
      };

}

#endif // HEMELB_IO_WRITERS_WRITER_H
