#ifndef HEME_IO_WRITER_H
#define HEME_IO_WRITER_H

#include "vis/colPixel.h"
#include "vis/colourPalette.h"

namespace heme 
{
  namespace io
  {
    class Writer
    {
    public:
    
      enum Separator {
	eol
      };
    
      // Overload << to write basic types and any necessary separators
      template <typename T>
	Writer& operator<< (T const & value) {
	_write(value);
	writeFieldSeparator();
	return *this;
      }
      // Special version for eol
      Writer& operator<< (enum Separator & value);
    
      void writePixel (vis::ColPixel *col_pixel_p,
		       vis::ColourPaletteFunction* colourPalette);
      //void (*colourPalette)(float value, float col[]));
    
      // Function to get the current position of writing in the stream.
      virtual unsigned int getCurrentStreamPosition() const = 0;
    
    protected:
    
      // Functions for formatting control
      virtual void writeFieldSeparator() = 0;
      virtual void writeRecordSeparator() = 0;
    
      // Methods to simply write (no separators) which are virtual and
      // hence must be overriden.
      virtual void _write(int const & value) = 0;
      virtual void _write(double const & value) = 0;
      virtual void _write(short const & value) = 0;
      virtual void _write(float const & value) = 0;
      virtual void _write(unsigned int const & value) = 0;
    
    };
  
  }
}

#endif // HEME_IO_WRITER_H
