// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_EXTRACTION_LOCALPROPERTYOUTPUT_H
#define HEMELB_EXTRACTION_LOCALPROPERTYOUTPUT_H

#include "extraction/IterableDataSource.h"
#include "extraction/PropertyOutputFile.h"
#include "io/TimePattern.h"
#include "lb/Lattices.h"
#include "net/mpi.h"
#include "net/IOCommunicator.h"
#include "net/MpiFile.h"

namespace hemelb::net { class IOCommunicator; }

namespace hemelb::extraction
{
    // Stores sufficient information to output property information
    // from this core.
    class LocalPropertyOutput
    {
    public:
      // Initialises a LocalPropertyOutput. Required so we can use
      // const reference types. Collective on the communicator.
      LocalPropertyOutput(std::shared_ptr<IterableDataSource> dataSource, const PropertyOutputFile& outputSpec,
			  net::IOCommunicator ioComms);

      // True if this property output should be written on the current iteration.
      bool ShouldWrite(unsigned long timestepNumber) const;

      // Returns the property output file object to be written.
      const PropertyOutputFile& GetOutputSpec() const;

      // Write this core's section of the data file. Only writes if
      // appropriate for the current iteration number.
      // Returns the path of the written file, iff writing occurred.
      std::optional<std::filesystem::path>
      Write(unsigned long timestepNumber, unsigned long totalSteps);

      inline std::string const& GetOffsetFileName() const {
          return offset_file_name;
      }

      // Write the offset file. Collective on the communicator.
      void WriteOffsetFile();

      // Returns the number of items written for the field.
      unsigned GetFieldLength(source::Type) const;

      auto OnIORank() const {
          return comms.OnIORank();
      }

    private:
      // How many sites does this MPI process write?
      std::uint64_t CountWrittenSitesOnRank();

      // How many bytes are written for a single site?
      std::uint64_t CalcSiteWriteLen(std::vector<OutputField> const& fields) const;

      // Make the XTR header
      std::vector<std::byte> PrepareHeader() const;

      // Open the file specified and write the header. Collective.
      void StartFile(std::string const& fn);

      // Our communicator
      net::IOCommunicator comms;

      // For single-timestep-per-file mode, hold the pattern we'll
      // pass to printf.
      io::TimePattern output_file_pattern;

      // The MPI file to write into.
      net::MpiFile outputFile;

      // The data source to use for file output.
      std::shared_ptr<IterableDataSource> dataSource;

      // PropertyOutputFile spec.
      PropertyOutputFile outputSpec;

      // How many local/global sites will be written
      std::uint64_t local_site_count;
      std::uint64_t global_site_count;

      // The length, in bytes, of the whole header
      std::uint64_t header_length;
      // The data that makes up the header (only used on rank 0)
      std::vector<std::byte> header_data;

      // The length, in bytes, of the local/global data write for one timestep
      std::uint64_t local_data_write_length;
      std::uint64_t global_data_write_length;

      // Where, in bytes, to begin writing into the file.
      std::uint64_t local_write_start;

      // Buffer to serialise into before writing to disk.
      std::vector<std::byte> buffer;

      // The file to write the offsets into.
      std::string offset_file_name;
    };
}

#endif // HEMELB_EXTRACTION_LOCALPROPERTYOUTPUT_H
