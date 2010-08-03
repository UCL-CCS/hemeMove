#ifndef __writer_h_
#define __writer_h_


#include "rt.h"

class Writer {
  
 public:
  // Functions to write basic types and any necessary separators
  template <typename T> void write(T const & value);
  
  void writePixel (ColPixel *col_pixel_p,
		   void (*colourPalette)(float value, float col[]));
  
  // Function to get the current position of writing in the stream.
  virtual unsigned int getCurrentStreamPosition() const = 0;
  
  // Functions for formatting control
  virtual void writeFieldSeparator() = 0;
  virtual void writeRecordSeparator() = 0;
  
 protected:
  // Constructor and member variable. This class should never be
  // instantiated directly, so the constructor is only available to
  // subclasses.
  
  Writer();
  
  // Methods to simply write (no separators) which are virtual and
  // hence must be overriden.
  virtual void _write(int const & value) = 0;
  virtual void _write(double const & value) = 0;
  virtual void _write(short const & value) = 0;
  virtual void _write(float const & value) = 0;
  virtual void _write(unsigned int const & value) = 0;
  
};

#endif //__writer_h_
