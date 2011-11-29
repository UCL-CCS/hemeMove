#ifndef HEMELB_REPORTING_PATHMANAGER_H
#define HEMELB_REPORTING_PATHMANAGER_H

#include <string>
#include "configuration/CommandLine.h"
#include "util/fileutils.h"
#include "log/Logger.h"
#include "configuration/SimConfig.h"
#include "io/writers/xdr/XdrFileWriter.h"
namespace hemelb
{
  namespace reporting
  {
    class PathManager
    {
      public:
        PathManager(configuration::CommandLine & commandLine, const bool &io, const int &processorCount);
        bool HasProblems() const
        {
          return(!ok);
        }
        const std::string & GetInputFile() const;
        const std::string & GetSnapshotDirectory() const;
        const std::string & GetImageDirectory() const;
        const std::string & GetReportPath() const;
        void SaveConfiguration(configuration::SimConfig * simConfig);
        void EmptyOutputDirectories();
        hemelb::io::writers::xdr::XdrFileWriter * XdrImageWriter(const long int time);
        const std::string SnapshotPath(unsigned long time) const;
      private:
        void GuessOutputDir();
        void InitialiseReport(const int & processorCount);
        std::string outputDir;
        std::string inputFile;
        std::string snapshotDirectory;
        std::string imageDirectory;
        std::string configLeafName;
        std::string timings_name;
        configuration::CommandLine &options;
        bool ok;
        bool doIo;
    };
  }
}

#endif //HEMELB_REPORTING_PATHMANAGER_H
