#include <fstream>

#include "io/AsciiFileWriter.h"

namespace hemelb
{
  namespace io
  {

    // Implement a constructor that opens the file and creates the Xdr
    // object to write to it.
    AsciiFileWriter::AsciiFileWriter(std::string fileName)
    {
      outStream = new std::ofstream(fileName.c_str());
    }

    // Destructor closes the file and cleans up the stream object.
    AsciiFileWriter::~AsciiFileWriter()
    {
      if (outStream != NULL)
      {
        // Since we know the member std::ostream* outStream is actually an
        // instance of std::ofstream, cast it to that so we can close it.
        std::ofstream* outFileStream = static_cast<std::ofstream *> (outStream);

        outFileStream->close();
        delete outFileStream;
      }
    }

  }
}
