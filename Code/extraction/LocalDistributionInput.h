// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_EXTRACTION_LOCALDISTRIBUTIONINPUT_H
#define HEMELB_EXTRACTION_LOCALDISTRIBUTIONINPUT_H

#include <optional>
#include <filesystem>

#include "extraction/IterableDataSource.h"
#include "extraction/InputField.h"
#include "io/readers/XdrMemReader.h"
#include "lb/Lattices.h"
#include "net/mpi.h"
#include "net/MpiFile.h"
#include "net/IOCommunicator.h"


namespace hemelb
{
  namespace net
  {
    class IOCommunicator;
  }
  namespace geometry
  {
    class FieldData;
  }
  namespace extraction
  {
    // Read this rank's part of a checkpoint file.
    class LocalDistributionInput
    {
      public:
      // Construct, but don't do any I/O.
      //
      // Take the path string by value since we will move it into a member
      // anyway.
      LocalDistributionInput(std::filesystem::path dataFilePath, std::optional<std::filesystem::path> maybeOffsetPath, const net::IOCommunicator& ioComms);

      // Open the file and load our part into the domain_type
      // instance.
      //
      // Time is optional, if not supplied will use the last one in
      // the file and will set the argument to that value.
      //
      // Requires the checkpoint have been saved with exactly the same
      // domain decomposition as currently running.
      void LoadDistribution(geometry::FieldData* latDat, std::optional<LatticeTimeStep>& initalTime);

    private:

      void ReadExtractionHeaders(net::MpiFile&, const unsigned NUMVECTORS);
      void ReadOffsets(const std::string&);

      const net::IOCommunicator& comms;

      // The path to the file to read from.
      std::filesystem::path filePath;
      std::filesystem::path offsetPath;

      InputField distField;
      uint64_t localStart;
      uint64_t localStop;
      uint64_t timestep;
      uint64_t allCoresWriteLength;
    };
  }
}

#endif // HEMELB_EXTRACTION_LOCALDISTRIBUTIONINPUT_H
