#include <iostream>

#include "io/AsciiStreamWriter.h"

namespace hemelb
{
  namespace io
  {

    // Method to get the current position in the stream.
    unsigned int AsciiStreamWriter::getCurrentStreamPosition() const
    {
      return (unsigned int) outStream->tellp();
    }

    // No field/record separators in XDR files
    void AsciiStreamWriter::writeFieldSeparator()
    {
      *outStream << ' ';
    }
    void AsciiStreamWriter::writeRecordSeparator()
    {
      *outStream << std::endl;
    }

    void AsciiStreamWriter::_write(int const & value)
    {
      this->_write<int> (value);
    }
    void AsciiStreamWriter::_write(double const & value)
    {
      this->_write<double> (value);
    }
    void AsciiStreamWriter::_write(float const & value)
    {
      this->_write<float> (value);
    }
    void AsciiStreamWriter::_write(short const & value)
    {
      this->_write<short> (value);
    }
    void AsciiStreamWriter::_write(unsigned int const & value)
    {
      this->_write<unsigned int> (value);
    }
    void AsciiStreamWriter::_write(long const & value)
    {
      this->_write<long> (value);
    }
  }
}
