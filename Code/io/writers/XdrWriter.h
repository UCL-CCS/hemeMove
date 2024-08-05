// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_IO_WRITERS_XDRWRITER_H
#define HEMELB_IO_WRITERS_XDRWRITER_H

#include "io/writers/Writer.h"
#include "io/XdrSerialisation.h"
#include "hassert.h"

namespace hemelb::io
{

	// Base XDR writer - mainly to give a common base class for
	// the generic ones.
	class XdrWriter : public Writer {
	protected:
	  void writeFieldSeparator() override;
	  void writeRecordSeparator() override;
	};

	// Main XDR writer class template
	//
	// ByteOutputIterator - the iterator we are going to use to
	// store bytes, should dereference to char (technically we
	// want:
	// 
	// *std::declval<ByteOutputIterator>() = std::declval<char>()
	// 
	// to be valid).
	//
	// Resource - a resource we might potentially use to store data
	template<std::output_iterator<std::byte> ByteOutputIterator,
		 typename Resource = xdr::Null>
	class XdrMetaWriter : public XdrWriter
	{
    protected:
	  // Let resource in, if need be
	  friend Resource;
	  Resource res;

	  ByteOutputIterator current;

	  // Some of the following could be an empty struct depending
	  // on the traits for the iterator
	  using boi_traits = xdr::boi_traits<ByteOutputIterator>;
	  typename boi_traits::counter_type bytes_written;
	  typename boi_traits::start_type start;
	  typename boi_traits::end_type end;

	  // Special constructor for difficult derived types to call.
	  //
	  // F and G are invokable to initialise the Resource and
	  // current BOI. Must have the following overloades:
	  //
	  // Resource F::operator()();
	  // ByteOperatorIterator G::operator()(Resource&);
	  //
	  // This is a bit of a hack to make Xdr{Vector, File}Writer
	  // easy...
	  template <typename F, typename G>
	  XdrMetaWriter(F f, G g) :
	    res(f()),
	    current(g(res)),
	    bytes_written(0),
	    start(current),
	    end()
	  {
	  }

	public:
	  // In the case of the Resource being Null or default
	  // constructible have two nice easy constructors.
	  //
	  // We have an iterator that supports indefinite insertion
	  explicit XdrMetaWriter(ByteOutputIterator start_) :
	    current(start_),
	    bytes_written(0),
	    start(start_),
	    end()
	  {
	  }
	  // We have a finite amount of space between start and end
	  XdrMetaWriter(ByteOutputIterator start_, ByteOutputIterator end_) :
	    current(start_),
	    bytes_written(0),
	    start(start_),
	    end(end_)
	  {
	    static_assert(xdr::is_rai<ByteOutputIterator>::value,
			  "Can't use this constructor for non-RandomAccessIterators!");
	  }

	  ~XdrMetaWriter() override  = default;

	  unsigned int getCurrentStreamPosition() const override {
	    return boi_traits::get_position(bytes_written, start, current);
	  }
	protected:
	  // Methods to simply write (no separators) which are virtual and
	  // hence must be overriden.
	  void _write(std::int16_t const& intToWrite) override {
	    write(intToWrite);
	  }
	  void _write(std::uint16_t const& uIntToWrite) override {
	    write(uIntToWrite);
	  }
	  void _write(std::int32_t const& intToWrite) override {
	    write(intToWrite);
	  }
	  void _write(std::uint32_t const& uIntToWrite) override {
	    write(uIntToWrite);
	  }
	  void _write(std::int64_t const& intToWrite) override {
	    write(intToWrite);
	  }
	  void _write(std::uint64_t const& uIntToWrite) override {
	    write(uIntToWrite);
	  }

	  void _write(double const& doubleToWrite) override {
	    write(doubleToWrite);
	  }
	  void _write(float const& floatToWrite) override {
	    write(floatToWrite);
	  }

	  void _write(const std::string& stringToWrite) override {
	    // The standard defines a string of n (numbered 0 through
	    // n-1) ASCII bytes to be the number n encoded as an
	    // unsigned integer (as described above), and followed by
	    // the n bytes of the string.  Byte m of the string always
	    // precedes byte m+1 of the string, and byte 0 of the
	    // string always follows the string's length.  If n is not
	    // a multiple of four, then the n bytes are followed by
	    // enough (0 to 3) residual zero bytes, r, to make the
	    // total byte count a multiple of four.
	    const uint32_t len = stringToWrite.size();
	    write(len);
        auto data = reinterpret_cast<std::byte const*>(stringToWrite.data());
        _write(std::span(data, data + len));
	  }

      // Like xdr_opaque
      void _write(std::span<const std::byte> bytes) override {
          const auto len = bytes.size();
          const auto req_nwords = (len - 1)/4 + 1;
          const auto space = req_nwords * 4;

          HASSERT(boi_traits::check_space(end, current, space));
          current = std::copy(bytes.begin(), bytes.end(), current);
          // Padding
          for (auto i = len; i != space; ++i)
              *current++ = std::byte{0};
          bytes_written += 4*req_nwords;
      }

	  template <typename T>
	  void write(T const& valToWrite) {
	    constexpr auto buf_size = xdr::xdr_serialised_size<T>();
	    // If we have an end then we must have space to store the serialised value
	    HASSERT(boi_traits::check_space(end, current, buf_size));

	    std::byte buf[buf_size];
	    xdr::xdr_serialise(valToWrite, buf);
	    current = std::copy(buf, buf + buf_size, current);
	    bytes_written += buf_size;
	  }
	};

    template <typename ItT>
    XdrMetaWriter<ItT> MakeXdrWriter(ItT start) {
      return XdrMetaWriter<ItT>(start);
    }
    template <typename ItT>
    XdrMetaWriter<ItT> MakeXdrWriter(ItT start, ItT end) {
      return XdrMetaWriter<ItT>(start, end);
    }

}

#endif //HEMELB_IO_WRITERS_XDRWRITER_H
