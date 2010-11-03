#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "constants.h"
#include "usage.h"
#include "fileutils.h"
#include "utilityFunctions.h"
#include "lb.h"

#include "SimConfig.h"
#include "SimulationMaster.h"

#include "steering/steering.h"

#include "vis/visthread.h"
#include "vis/ColourPalette.h"

#include "debug/Debugger.h"

FILE *timings_ptr;

int main(int argc, char *argv[])
{
  // main function needed to perform the entire simulation. Some
  // simulation paramenters and performance statistics are outputted on
  // standard output

  SimulationMaster lMaster = SimulationMaster(argc, argv);

  // This is currently where all default command-line arguments are.
  std::string lInputFile = "config.xml";
  std::string lOutputDir = "";
  unsigned int lSnapshotsPerCycle = 10;
  unsigned int lImagesPerCycle = 10;
  unsigned int lSteeringSessionId = 1;

  // There should be an odd number of arguments since the parameters occur in pairs.
  if ( (argc % 2) == 0)
  {
    if (lMaster.GetNet()->IsCurrentProcTheIOProc())
    {
      Usage::printUsage(argv[0]);
    }
    lMaster.GetNet()->Abort();
  }

  // All arguments are parsed in pairs, one is a "-<paramName>" type, and one
  // is the <parametervalue>.
  for (int ii = 1; ii < argc; ii += 2)
  {
    char* lParamName = argv[ii];
    char* lParamValue = argv[ii + 1];
    if (strcmp(lParamName, "-in") == 0)
    {
      lInputFile = std::string(lParamValue);
    }
    else if (strcmp(lParamName, "-out") == 0)
    {
      lOutputDir = std::string(lParamValue);
    }
    else if (strcmp(lParamName, "-s") == 0)
    {
      char * dummy;
      lSnapshotsPerCycle = (unsigned int) (strtoul(lParamValue, &dummy, 10));
    }
    else if (strcmp(lParamName, "-i") == 0)
    {
      char *dummy;
      lImagesPerCycle = (unsigned int) (strtoul(lParamValue, &dummy, 10));
    }
    else if (strcmp(lParamName, "-ss") == 0)
    {
      char *dummy;
      lSteeringSessionId = (unsigned int) (strtoul(lParamValue, &dummy, 10));
    }
    else
    {
      if (lMaster.GetNet()->IsCurrentProcTheIOProc())
      {
        Usage::printUsage(argv[0]);
      }
      lMaster.GetNet()->Abort();
    }

  }

  hemelb::SimConfig *lSimulationConfig =
      hemelb::SimConfig::Load(lInputFile.c_str());

  unsigned long lLastForwardSlash = lInputFile.rfind('/');
  if (lOutputDir.length() == 0)
  {
    lOutputDir = ( (lLastForwardSlash == std::string::npos)
      ? "./"
      : lInputFile.substr(0, lLastForwardSlash)) + "results";
  }

  double total_time = hemelb::util::myClock();
  char snapshot_directory[256];
  char image_directory[256];
  // Create directory path for the output images
  strcpy(image_directory, lOutputDir.c_str());
  strcat(image_directory, "/Images/");
  //Create directory path for the output snapshots
  strcpy(snapshot_directory, lOutputDir.c_str());
  strcat(snapshot_directory, "/Snapshots/");
  // Actually create the directories.

  if (lMaster.GetNet()->IsCurrentProcTheIOProc())
  {
    if (hemelb::util::DoesDirectoryExist(lOutputDir.c_str()))
    {
      printf("\nOutput directory \"%s\" already exists. Exiting.\n\n",
             lOutputDir.c_str());
      lMaster.GetNet()->Abort();
    }
    hemelb::util::MakeDirAllRXW(lOutputDir.c_str());
    hemelb::util::MakeDirAllRXW(image_directory);
    hemelb::util::MakeDirAllRXW(snapshot_directory);
    // Save the computed config out to disk in the output directory so we have
    // a record of the total state used.
    std::string lFileNameComponent = std::string( (lLastForwardSlash
        == std::string::npos)
      ? lInputFile
      : lInputFile.substr(lLastForwardSlash));
    char *lConfigCopyName = new char[lOutputDir.length()
        + lFileNameComponent.length() + 3];
    sprintf(lConfigCopyName, "%s/%s", lOutputDir.c_str(),
            lFileNameComponent.c_str());
    lSimulationConfig->Save(lConfigCopyName);
    delete lConfigCopyName;
    char timings_name[256];
    char procs_string[256];
    sprintf(procs_string, "%i", lMaster.GetNet()->mProcessorCount);
    strcpy(timings_name, lOutputDir.c_str());
    strcat(timings_name, "/timings");
    strcat(timings_name, procs_string);
    strcat(timings_name, ".asc");
    timings_ptr = fopen(timings_name, "w");
    fprintf(timings_ptr,
            "***********************************************************\n");
    fprintf(timings_ptr, "Opening config file:\n %s\n", lInputFile.c_str());
  }

  lMaster.Initialise(lSimulationConfig, (int) lSteeringSessionId, timings_ptr);

  lMaster.RunSimulation(timings_ptr, lSimulationConfig, total_time,
                        image_directory, snapshot_directory,
                        lSnapshotsPerCycle, lImagesPerCycle);

  if (lMaster.GetNet()->IsCurrentProcTheIOProc())
  {
    fclose(timings_ptr);
  }

  delete lSimulationConfig;
  delete hemelb::vis::controller;

  return (0);
}

