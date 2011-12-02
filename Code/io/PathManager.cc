#include "PathManager.h"
#include <sstream>
namespace hemelb
{
  namespace io
  {
    PathManager::PathManager(configuration::CommandLine &commandLine,
                             const bool & io,
                             const int & processorCount) :
        options(commandLine), ok(false), doIo(io)
    {

      inputFile = options.GetInputFile();
      outputDir = options.GetOutputDir();

      GuessOutputDir();

      imageDirectory = outputDir + "/Images/";
      snapshotDirectory = outputDir + "/Snapshots/";

      if (doIo)
      {
        if (hemelb::util::DoesDirectoryExist(outputDir.c_str()))
        {
          hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("\nOutput directory \"%s\" already exists. Exiting.",
                                                                              outputDir.c_str());
          return;
        }

        hemelb::util::MakeDirAllRXW(outputDir);
        hemelb::util::MakeDirAllRXW(imageDirectory);
        hemelb::util::MakeDirAllRXW(snapshotDirectory);
        std::stringstream reportNameStream;
        reportNameStream << outputDir << "/timings" << processorCount << ".asc" << std::flush;
        reportName = reportNameStream.str();
      }

      ok = true;
    }

    const std::string & PathManager::GetInputFile() const
    {
      return inputFile;
    }
    const std::string & PathManager::GetSnapshotDirectory() const
    {
      return snapshotDirectory;
    }
    const std::string & PathManager::GetImageDirectory() const
    {
      return imageDirectory;
    }
    const std::string & PathManager::GetReportPath() const
    {
      return reportName;
    }

    void PathManager::EmptyOutputDirectories()
    {
      hemelb::util::DeleteDirContents(snapshotDirectory);
      hemelb::util::DeleteDirContents(imageDirectory);
    }

    hemelb::io::writers::xdr::XdrFileWriter * PathManager::XdrImageWriter(const long int time)
    {
      char filename[255];
      snprintf(filename, 255, "%08li.dat", time);
      return (new hemelb::io::writers::xdr::XdrFileWriter(imageDirectory + std::string(filename)));
    }

    const std::string PathManager::SnapshotPath(unsigned long time) const
    {
      char snapshotFilename[255];
      snprintf(snapshotFilename, 255, "snapshot_%06li.dat", time);
      return (snapshotDirectory + std::string(snapshotFilename)); // by copy
    }

    void PathManager::SaveConfiguration(configuration::SimConfig * simConfig)
    {
      if (doIo)
      {
        simConfig->Save(outputDir + "/" + configLeafName);
      }
    }

    void PathManager::GuessOutputDir()
    {
      unsigned long lLastForwardSlash = inputFile.rfind('/');
      if (lLastForwardSlash == std::string::npos)
      {
        // input file supplied is in current folder
        configLeafName = inputFile;
        if (outputDir.length() == 0)
        {
          // no output dir given, defaulting to local.
          outputDir = "./results";
        }
      }
      else
      {
        // input file supplied is a path to the input file
        configLeafName = inputFile.substr(lLastForwardSlash);
        if (outputDir.length() == 0)
        {
          // no output dir given, defaulting to location of input file.
          // note substr is end-exclusive and start-inclusive
          outputDir = inputFile.substr(0, lLastForwardSlash + 1) + "results";
        }
      }
    }
  }
}

