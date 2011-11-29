#ifndef HEMELB_IO_XDRFILEWRITER_H
#define HEMELB_IO_XDRFILEWRITER_H

#include <stdio.h>
#include <string>

#include "io/writers/xdr/XdrWriter.h"

namespace hemelb
{
  namespace io
  {
    namespace writers
    {
      namespace xdr
      {

        // Class to write Xdr to a file. The actual write functions are implemented in the base class, XdrWriter.
        class XdrFileWriter : public XdrWriter
        {

            // Implement the constructor and destructor to deal with the FILE
            // and XDR objects.
          public:
            XdrFileWriter(const std::string fileName, const std::string mode = "w");
            ~XdrFileWriter();

          private:
            FILE *myFile;

        };

      } // namespace xdr
    } // namespace writers
  }
}
#endif // HEMELB_IO_XDRFILEWRITER_H
