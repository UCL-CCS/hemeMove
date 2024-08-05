// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_IO_XDRSERIALISATION_H
#define HEMELB_IO_XDRSERIALISATION_H

#include <type_traits>
#include <iterator>
#include <optional>

#include <arpa/inet.h>

namespace hemelb::io::xdr {
	  // Really implementation details of the XDR coding

	  // Return the XDR representation size of a type
	  template <typename T>
	  constexpr size_t xdr_serialised_size() {
	    return 4 * ((sizeof(T)-1)/4 + 1);
	  }

	  // Overload set for XDR primitives to serialise a value into a buffer
	  // Basic signature is :
	  //  void xdr_serialise(const T& value, char* buffer);

	  // Signed and unsigned integers of 32 bits or less
	  template<std::integral T>
      requires (sizeof(T) <= 4)
	  void xdr_serialise(const T& val, std::byte* dest_buf)
	  {
	    auto out = reinterpret_cast<uint32_t*>(dest_buf);
	    out[0] = htonl(val);
	  }
	    
	  // Signed and unsigned 64 bit ints
	  template<std::integral T>
	  requires (sizeof(T) == 8)
	  void xdr_serialise(const T& val, std::byte* dest_buf)
	  {
	    auto out = reinterpret_cast<uint32_t*>(dest_buf);
	    auto msw = uint32_t(val >> 32);
	    auto lsw = uint32_t(val);
	    out[0] = htonl(msw);
	    out[1] = htonl(lsw);
	  }

	  // 32 bit floats
	  inline void xdr_serialise(const float& val, std::byte* dest_buf)
	  {
	    auto data = reinterpret_cast<const uint32_t*>(&val);
	    auto out = reinterpret_cast<uint32_t*>(dest_buf);
	    out[0] = htonl(*data);
	  }
	  // 64 bit doubles
	  inline void xdr_serialise(const double& val, std::byte* dest_buf)
	  {
	    auto data = reinterpret_cast<const uint64_t*>(&val);
	    xdr_serialise(*data, dest_buf);
	  }
	  // End of xdr_serialise overloads

	  // Overload set for XDR primitives to deserialise a value from a buffer
	  // Basic signature is :
	  //  void xdr_deserialise(T& dest, const char* buffer);
	  // Signed and unsigned integers of 32 bits or less
	  template<std::integral T>
	  requires (sizeof(T) <= 4)
	  void xdr_deserialise(T& val, const std::byte* src_buf)
	  {
	    auto in = reinterpret_cast<const uint32_t*>(src_buf);
	    val = ntohl(*in);
	  }
	    
	  // Signed and unsigned 64 bit ints
	  template<std::integral T>
	  requires (sizeof(T) == 8)
	  void xdr_deserialise(T& val, const std::byte* src_buf)
	  {
	    auto in = reinterpret_cast<const uint32_t*>(src_buf);
	    uint32_t msw = ntohl(in[0]);
	    uint32_t lsw = ntohl(in[1]);
	    val = T(msw) << 32 ^ lsw;
	  }

	  // 32 bit floats
	  inline void xdr_deserialise(float& val, const std::byte* src_buf)
	  {
	    auto in = reinterpret_cast<const uint32_t*>(src_buf);
	    auto out = reinterpret_cast<uint32_t*>(&val);
	    *out = ntohl(*in);
	  }

	  // 64 bit doubles
	  inline void xdr_deserialise(double& val, const std::byte* src_buf)
	  {
	    auto out = reinterpret_cast<uint64_t*>(&val);
	    xdr_deserialise(*out, src_buf);
	  }
	  // End of xdr_serialise overloads
	  
	  
	  // A completely empty struct to be used as a placeholder
	  struct Null {
	    template <typename... Ts>
	    constexpr Null(Ts...) {
	    }
	    inline Null& operator+=(size_t) { return *this; }
	    constexpr bool operator!() const {
	      return true;
	    }
	    constexpr char* operator*() const {
	      return nullptr;
	    }
	  };

	  // Is the parameter a random access iterator?
	  // Answer in static member value
	  template <typename I>
	  struct is_rai {
	    static constexpr bool value = std::is_same<
	      typename std::iterator_traits<I>::iterator_category,
	      std::random_access_iterator_tag
	      >::value;
	  };

	  template <typename ItT>
	  struct boi_traits_nonaddable
	  {
	    static constexpr bool can_add = false;
	    // No need to track start
	    using start_type = Null;
	    // Or end - indeed without `operator<=` it won't work
	    using end_type = Null;
	    // Instead we count the number of bytes we write
	    using counter_type = size_t;
	    // Assume that we can write forever - e.g. to a file
	    static constexpr bool check_space(const end_type& end, const ItT& current, size_t buf_size) {
	      return true;
	    }
	    // The current position is just the number of bytes written
	    static unsigned get_position(const counter_type& bytes_written,
					 const start_type& start,
					 const ItT& current) {
	      return bytes_written;
	    }
	  };

	  // Here we specialise for the case when we do have random
	  // access iterators
	  template <typename ItT>
	  struct boi_traits_addable
	  {
	    static constexpr bool can_add = true;
	    // Start is just the same (fingers crossed iterators
	    // aren't invalidated...)
	    using start_type = ItT;
	    // End is an optional iter
	    using end_type = std::optional<ItT>;

	    // No counter needed (Null supports operator+=(size_t))
	    using counter_type = Null;
	    // Check that we have (unlimited) space
	    static bool check_space(const end_type& end, const ItT& current, size_t buf_size) {
	      return !end || (current + buf_size <= *end);
	    }
	    // The current position is how far we got from the start
	    static unsigned get_position(const counter_type& bytes_written,
					 const start_type& start,
					 const ItT& current) {
	      return current - start;
	    }
	  };

	  // Byte Output Iterator traits - to support XdrMetaWriter
	  //
	  // We have two cases: when the iterators support random
	  // access and when they don't.
	  template <typename I>
	  using boi_traits = typename std::conditional<
	    is_rai<I>::value,
	    boi_traits_addable<I>,
	    boi_traits_nonaddable<I>
	    >::type;

}

#endif
