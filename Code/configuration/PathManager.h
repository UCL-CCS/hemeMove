// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_CONFIGURATION_PATHMANAGER_H
#define HEMELB_CONFIGURATION_PATHMANAGER_H

#include <filesystem>
#include <string>
#include "log/Logger.h"

namespace hemelb::configuration
{
    class CommandLine;
    class SimConfig;

    /**
     * Manage the input and output file system locations for HemeLB reports, extracted data, and input xml.
     */
    class PathManager
    {
      using path = std::filesystem::path;
      public:
        /**
         * During construction, the path manager will guess an appropriate path from the input file path given at the command line,
         * if one was not specified explicitly.
         * @param commandLine The command-line parser object instantiated from the arguments to the C main function.
         * @param io Set to true if this MPI node is the I/O process and files should be written from this node.
         * @param processorCount The total count of processors, used in generating file-names for report files.
         */
        PathManager(const configuration::CommandLine & commandLine, const bool &io,
                    const int &processorCount);

        /**
         * A local or full path to the input xml configuration file.
         * @return A local or full path to the input xml configuration file.
         */
        [[nodiscard]] const path& GetInputFile() const;
        /**
         * Gets the path to the file where colloid output should be written
         * @return
         */
        [[nodiscard]] const path& GetColloidPath() const;
        /**
         * Path to where a run report file should be created.
         * @return Reference to path to where a run report file should be created.
         */
        [[nodiscard]] const path& GetReportPath() const;

        /**
         * Return the path that property extraction output should go to.
         * @return
         */
        [[nodiscard]] const path& GetDataExtractionPath() const;

        [[nodiscard]] const path& GetCheckpointPath() const;

        /**
         * Create a subdirectory inside the RBC output directory and return its path
         * @param subdirectoryName Name of the subdirectory to be created
         * @return Path to the newly created subdirectory
         */
        [[nodiscard]] path GetRBCOutputPathWithSubdir(std::string const& subdirectoryName, bool create, bool allow_existing) const;

      private:
        path outputDir;
        path inputFile;
        path colloidFile;
        path extractionDir;
        path cpDir;
        const configuration::CommandLine &options;
        bool doIo; //! Am I the input/output node?
        path rbcDir; //! Path for RBC output
    };

}

#endif //HEMELB_CONFIGURATION_PATHMANAGER_H
