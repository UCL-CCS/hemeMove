#include "SimulationMaster.h"
#include "SimConfig.h"

#include "io/XdrFileWriter.h"
#include "util/utilityFunctions.h"
#include "geometry/LatticeData.h"
#include "debug/Debugger.h"
#include "util/fileutils.h"
#include "log/Logger.h"
#include "lb/HFunction.h"

#include "topology/NetworkTopology.h"

#include <map>
#include <limits>
#include <stdlib.h>

/**
 * Constructor for the SimulationMaster class
 *
 * Initialises member variables including the network topology
 * object.
 */
SimulationMaster::SimulationMaster(int iArgCount, char *iArgList[])
{
  // Initialise the network discovery. If this fails, abort.
  bool lTopologySuccess = true;
  hemelb::topology::NetworkTopology::Instance()->Init(&iArgCount, &iArgList, &lTopologySuccess);

  if (!lTopologySuccess)
  {
    hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::OnePerCore>("Couldn't get machine information for this network topology. Aborting.\n");
    Abort();
  }

  mCreationTime = hemelb::util::myClock();

  hemelb::debug::Debugger::Init(iArgList[0]);

  mLatDat = NULL;

  mLbm = NULL;
  steeringCpt = NULL;
  mVisControl = NULL;
  mSimulationState = NULL;

  // This is currently where all default command-line arguments are.
  inputFile = "config.xml";
  outputDir = "";
  snapshotsPerCycle = 10;
  imagesPerCycle = 10;
  steeringSessionId = 1;

  mMPISendTime = 0.0;
  mMPIWaitTime = 0.0;
  mSnapshotTime = 0.0;

  mImagesWritten = 0;
  mSnapshotsWritten = 0;
  ParseArguments(iArgCount,iArgList);
  simConfig = hemelb::SimConfig::Load(inputFile.c_str());
  SetupReporting();
  Initialise();

}

void SimulationMaster::ParseArguments(int argc, char **argv){

  // There should be an odd number of arguments since the parameters occur in pairs.
    if ( (argc % 2) == 0)
    {
      if (IsCurrentProcTheIOProc())
      {
        PrintUsage();
      }
      Abort();
    }

  // All arguments are parsed in pairs, one is a "-<paramName>" type, and one
  // is the <parametervalue>.
  for (int ii = 1; ii < argc; ii += 2)
  {
    char* lParamName = argv[ii];
    char* lParamValue = argv[ii + 1];
    if (strcmp(lParamName, "-in") == 0)
    {
      inputFile = std::string(lParamValue);
    }
    else if (strcmp(lParamName, "-out") == 0)
    {
      outputDir = std::string(lParamValue);
    }
    else if (strcmp(lParamName, "-s") == 0)
    {
      char * dummy;
      snapshotsPerCycle = (unsigned int) (strtoul(lParamValue, &dummy, 10));
    }
    else if (strcmp(lParamName, "-i") == 0)
    {
      char *dummy;
      imagesPerCycle = (unsigned int) (strtoul(lParamValue, &dummy, 10));
    }
    else if (strcmp(lParamName, "-ss") == 0)
    {
      char *dummy;
      steeringSessionId = (unsigned int) (strtoul(lParamValue, &dummy, 10));
    }
    else
    {
      if (IsCurrentProcTheIOProc())
      {
        PrintUsage();
      }
      Abort();
    }
  }
}

void SimulationMaster::PrintUsage()
{
  printf("-!-!-!-!-!-!-!-!-!-!-!-!");
  printf("Correct usage: hemelb [-<Parameter Name> <Parameter Value>]* \n");
  printf("Parameter name and significance:\n");
  printf("-in \t Path to the configuration xml file (default is config.xml)\n");
  printf("-out \t Path to the output folder (default is based on input file, e.g. config_xml_results)\n");
  printf("-s \t Number of snapshots to take per cycle (default 10)\n");
  printf("-i \t Number of images to create per cycle (default is 10)\n");
  printf("-ss \t Steering session identifier (default is 1)\n");
  printf("-!-!-!-!-!-!-!-!-!-!-!-!");
}

void SimulationMaster::SetupReporting()
{

  std::string configLeafName;
  unsigned long lLastForwardSlash = inputFile.rfind('/');
  if (lLastForwardSlash == std::string::npos) {
    // input file supplied is in current folder
    configLeafName= inputFile;
    if (outputDir.length() == 0) {
      // no output dir given, defaulting to local.
      outputDir="./results";
    }
  } else {
    // input file supplied is a path to the input file
    configLeafName=  inputFile.substr(lLastForwardSlash);
    if (outputDir.length() == 0) {
      // no output dir given, defaulting to location of input file.
     outputDir=inputFile.substr(0, lLastForwardSlash)+"results";
    }
  }

  imageDirectory = outputDir + "/Images/";
  snapshotDirectory = outputDir + "/Snapshots/";

  if (IsCurrentProcTheIOProc())
    {
      if (hemelb::util::DoesDirectoryExist(outputDir.c_str()))
      {
        hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("\nOutput directory \"%s\" already exists. Exiting.",
                                                                            outputDir.c_str());
        Abort();
        return;
      }


      hemelb::util::MakeDirAllRXW(outputDir);
      hemelb::util::MakeDirAllRXW(imageDirectory);
      hemelb::util::MakeDirAllRXW(snapshotDirectory);

      simConfig->Save(outputDir + "/" + configLeafName);

      char timings_name[256];
      char procs_string[256];

      sprintf(procs_string, "%i", GetProcessorCount());
      strcpy(timings_name, outputDir.c_str());
      strcat(timings_name, "/timings");
      strcat(timings_name, procs_string);
      strcat(timings_name, ".asc");

      mTimingsFile = fopen(timings_name, "w");
      fprintf(mTimingsFile, "***********************************************************\n");
      fprintf(mTimingsFile, "Opening config file:\n %s\n", inputFile.c_str());
    }
}

/**
 * Destructor for the SimulationMaster class.
 *
 * Deallocates dynamically allocated memory to contained classes.
 */
SimulationMaster::~SimulationMaster()
{

  if (IsCurrentProcTheIOProc())
  {
   fclose(mTimingsFile);
  }
  if (hemelb::topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
  {
    delete imageSendCpt;
  }

  if (mLatDat != NULL)
  {
    delete mLatDat;
  }

  if (mLbm != NULL)
  {
    delete mLbm;
  }

  if (mInletValues != NULL)
  {
    delete mInletValues;
  }

  if (mOutletValues != NULL)
  {
    delete mOutletValues;
  }

  if (network != NULL)
  {
    delete network;
  }
  if (steeringCpt != NULL)
  {
    delete steeringCpt;
  }

  if (mVisControl != NULL)
  {
    delete mVisControl;
  }

  if (mStabilityTester != NULL)
  {
    delete mStabilityTester;
  }

  if (mEntropyTester != NULL)
  {
    delete mEntropyTester;
  }

  if (mSimulationState != NULL)
  {
    delete mSimulationState;
  }

  if (mUnits != NULL)
  {
    delete mUnits;
  }
  delete simConfig;
}

/**
 * Returns true if the current processor is the dedicated I/O
 * processor.
 */
bool SimulationMaster::IsCurrentProcTheIOProc()
{
  return hemelb::topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc();
}

/**
 * Returns the number of processors involved in the simulation.
 */
int SimulationMaster::GetProcessorCount()
{
  return hemelb::topology::NetworkTopology::Instance()->GetProcessorCount();
}

/**
 * Initialises various elements of the simulation if necessary - steering,
 * domain decomposition, LBM and visualisation.
 */
void SimulationMaster::Initialise()
{

  hemelb::log::Logger::Log<hemelb::log::Warning, hemelb::log::Singleton>("Beginning Initialisation.");

  mSimulationState = new hemelb::lb::SimulationState(simConfig->StepsPerCycle,
                                                     simConfig->NumCycles);

  hemelb::site_t mins[3], maxes[3];
  // TODO The way we initialise LbmParameters is not great.
  hemelb::lb::LbmParameters params(1000, 0.1);
  hemelb::site_t totalFluidSites;

  hemelb::log::Logger::Log<hemelb::log::Warning, hemelb::log::Singleton>("Initialising LatticeData.");
  mLatDat
      = new hemelb::geometry::LatticeData(hemelb::steering::SteeringComponent::RequiresSeparateSteeringCore(),
                                          &totalFluidSites,
                                          mins,
                                          maxes,
                                          hemelb::topology::NetworkTopology::Instance()->FluidSitesOnEachProcessor,
                                          &params,
                                          simConfig,
                                          &mFileReadTime,
                                          &mDomainDecompTime);

  hemelb::log::Logger::Log<hemelb::log::Warning, hemelb::log::Singleton>("Initialising LBM.");
  mLbm = new hemelb::lb::LBM(simConfig, &mNet, mLatDat, mSimulationState);
  mLbm->SetSiteMinima(mins);
  mLbm->SetSiteMaxima(maxes);

  mLbm->GetLbmParams()->StressType = params.StressType;
  mLbm->SetTotalFluidSiteCount(totalFluidSites);

  // Initialise and begin the steering.
  if (hemelb::topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
  {
    network = new hemelb::steering::Network(steeringSessionId);
  }
  else
  {
    network = NULL;
  }

  mStabilityTester = new hemelb::lb::StabilityTester(mLatDat, &mNet, mSimulationState);

  if (hemelb::log::Logger::ShouldDisplay<hemelb::log::Debug>())
  {
    int typesTested[1] = { 0 };
    mEntropyTester
        = new hemelb::lb::EntropyTester(typesTested, 1, mLatDat, &mNet, mSimulationState);
  }
  else
  {
    mEntropyTester = NULL;
  }

  double seconds = hemelb::util::myClock();
  hemelb::site_t* lReceiveTranslator = mNet.Initialise(mLatDat);
  mNetInitialiseTime = hemelb::util::myClock() - seconds;

  hemelb::log::Logger::Log<hemelb::log::Warning, hemelb::log::Singleton>("Initialising visualisation controller.");
  mVisControl = new hemelb::vis::Control(mLbm->GetLbmParams()->StressType,
                                         &mNet,
                                         mSimulationState,
                                         mLatDat);

  if (hemelb::topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
  {
    imageSendCpt = new hemelb::steering::ImageSendComponent(mLbm,
                                                            mSimulationState,
                                                            mVisControl,
                                                            mLbm->GetLbmParams(),
                                                            network);
  }

  mUnits = new hemelb::util::UnitConverter(mLbm->GetLbmParams(),
                                           mSimulationState,
                                           mLatDat->GetVoxelSize());

  mInletValues
      = new hemelb::lb::boundaries::BoundaryValues(hemelb::geometry::LatticeData::INLET_TYPE,
                                                   mLatDat,
                                                   simConfig->Inlets,
                                                   mSimulationState,
                                                   mUnits);

  mOutletValues
      = new hemelb::lb::boundaries::BoundaryValues(hemelb::geometry::LatticeData::OUTLET_TYPE,
                                                   mLatDat,
                                                   simConfig->Outlets,
                                                   mSimulationState,
                                                   mUnits);

  mLbm->Initialise(lReceiveTranslator, mVisControl, mInletValues, mOutletValues, mUnits);

  steeringCpt = new hemelb::steering::SteeringComponent(network,
                                                        mVisControl,
                                                        mLbm,
                                                        &mNet,
                                                        mSimulationState,
                                                        simConfig,
                                                        mUnits);

  // Read in the visualisation parameters.
  mLbm->ReadVisParameters();
}

/**
 * Begin the simulation.
 */
void SimulationMaster::RunSimulation()
{
  hemelb::log::Logger::Log<hemelb::log::Warning, hemelb::log::Singleton>("Beginning to run simulation.");

  double simulation_time = hemelb::util::myClock();
  bool is_unstable = false;
  int total_time_steps = 0;

  // TODO ugh.
  unsigned int
      snapshots_period =
          (snapshotsPerCycle == 0)
            ? 1000000000
            : hemelb::util::NumericalFunctions::max(1U,
                                                    (unsigned int) (mSimulationState->GetTimeStepsPerCycle()
                                                        / snapshotsPerCycle));

  unsigned int
      images_period =
          (imagesPerCycle == 0)
            ? 1000000000
            : hemelb::util::NumericalFunctions::max(1U,
                                                    (unsigned int) (mSimulationState->GetTimeStepsPerCycle()
                                                        / imagesPerCycle));

  bool is_finished = false;
  hemelb::lb::Stability stability = hemelb::lb::Stable;

  typedef std::multimap<unsigned long, unsigned long> mapType;

  mapType snapshotsCompleted;
  mapType networkImagesCompleted;

  std::vector<hemelb::net::IteratedAction*> actors;
  actors.push_back(mLbm);
  actors.push_back(mInletValues);
  actors.push_back(mOutletValues);
  actors.push_back(steeringCpt);
  actors.push_back(mStabilityTester);
  if (mEntropyTester != NULL)
  {
    actors.push_back(mEntropyTester);
  }
  actors.push_back(mVisControl);

  if (hemelb::topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
  {
    actors.push_back(network);
  }

  for (; mSimulationState->GetTimeStepsPassed() <= mSimulationState->GetTotalTimeSteps()
      && !is_finished; mSimulationState->Increment())
  {
    ++total_time_steps;

    bool write_snapshot_image = ( (mSimulationState->GetTimeStep() % images_period) == 0)
      ? true
      : false;

    // Make sure we're rendering if we're writing this iteration.
    if (write_snapshot_image)
    {
      snapshotsCompleted.insert(std::pair<unsigned long, unsigned long>(mVisControl->Start(),
                                                                        mSimulationState->GetTimeStepsPassed()));
    }

    if (mSimulationState->GetDoRendering())
    {
      networkImagesCompleted.insert(std::pair<unsigned long, unsigned long>(mVisControl->Start(),
                                                                            mSimulationState->GetTimeStepsPassed()));
      mSimulationState->SetDoRendering(false);
    }

    /* In the following two if blocks we do the core magic to ensure we only Render
     when (1) we are not sending a frame or (2) we need to output to disk */

    /* for debugging purposes we want to ensure we capture the variables in a single
     instant of time since variables might be altered by the thread half way through?
     This is to be done. */

    bool render_for_network_stream = false;
    if (hemelb::topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
    {
      render_for_network_stream = imageSendCpt->ShouldRenderNewNetworkImage();
      steeringCpt->readyForNextImage = render_for_network_stream;
    }

    if (mSimulationState->GetTimeStep() % 100 == 0)
    {
      hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("time step %i render_network_stream %i write_snapshot_image %i rendering %i",
                                                                          mSimulationState->GetTimeStep(),
                                                                          render_for_network_stream,
                                                                          write_snapshot_image,
                                                                          mSimulationState->GetDoRendering());
    }

    // Cycle.

    {
      for (std::vector<hemelb::net::IteratedAction*>::iterator it = actors.begin(); it
          != actors.end(); ++it)
      {
        (*it)->RequestComms();
      }

      mNet.Receive();

      {
        for (std::vector<hemelb::net::IteratedAction*>::iterator it = actors.begin(); it
            != actors.end(); ++it)
        {
          (*it)->PreSend();
        }
      }

      {
        double lPreSendTime = hemelb::util::myClock();
        mNet.Send();
        mMPISendTime += (hemelb::util::myClock() - lPreSendTime);
      }

      {
        for (std::vector<hemelb::net::IteratedAction*>::iterator it = actors.begin(); it
            != actors.end(); ++it)
        {
          (*it)->PreReceive();
        }
      }

      {
        double lPreWaitTime = hemelb::util::myClock();
        mNet.Wait();
        mMPIWaitTime += (hemelb::util::myClock() - lPreWaitTime);
      }

      {
        for (std::vector<hemelb::net::IteratedAction*>::iterator it = actors.begin(); it
            != actors.end(); ++it)
        {
          (*it)->PostReceive();
        }
      }

      for (std::vector<hemelb::net::IteratedAction*>::iterator it = actors.begin(); it
          != actors.end(); ++it)
      {
        (*it)->EndIteration();
      }
    }

    stability = hemelb::lb::Stable;

    if (mSimulationState->GetStability() == hemelb::lb::Unstable)
    {
      hemelb::util::DeleteDirContents(snapshotDirectory);
      hemelb::util::DeleteDirContents(imageDirectory);

      for (std::vector<hemelb::net::IteratedAction*>::iterator it = actors.begin(); it
          != actors.end(); ++it)
      {
        (*it)->Reset();
      }

#ifndef NO_STREAKLINES
      mVisControl->Reset();
#endif

      hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("restarting: period: %i\n",
                                                                          mSimulationState->GetTimeStepsPerCycle());

      snapshots_period
          = (snapshotsPerCycle == 0)
            ? 1000000000
            : hemelb::util::NumericalFunctions::max(1U,
                                                    (unsigned int) (mSimulationState->GetTimeStepsPerCycle()
                                                        / snapshotsPerCycle));

      images_period
          = (imagesPerCycle == 0)
            ? 1000000000
            : hemelb::util::NumericalFunctions::max(1U,
                                                    (unsigned int) (mSimulationState->GetTimeStepsPerCycle()
                                                        / imagesPerCycle));

      mSimulationState->Reset();
      continue;
    }
    mLbm->UpdateInletVelocities(mSimulationState->GetTimeStep());

#ifndef NO_STREAKLINES
    mVisControl->ProgressStreaklines(mSimulationState->GetTimeStep(),
                                     mSimulationState->GetTimeStepsPerCycle());
#endif

    if (snapshotsCompleted.count(mSimulationState->GetTimeStepsPassed()) > 0)
    {
      for (std::multimap<unsigned long, unsigned long>::const_iterator it =
          snapshotsCompleted.find(mSimulationState->GetTimeStepsPassed()); it
          != snapshotsCompleted.end() && it->first == mSimulationState->GetTimeStepsPassed(); ++it)
      {
        mImagesWritten++;

        if (hemelb::topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
        {
          char image_filename[255];
          snprintf(image_filename, 255, "%08li.dat", 1 + ( (it->second - 1)
              % mSimulationState->GetTimeStepsPerCycle()));
          hemelb::io::XdrFileWriter writer = hemelb::io::XdrFileWriter(imageDirectory
              + std::string(image_filename));

          const hemelb::vis::PixelSet<hemelb::vis::ResultPixel>* result =
              mVisControl->GetResult(it->second);

          mVisControl ->WriteImage(&writer,
                                   *result,
                                   &mVisControl->mDomainStats,
                                   &mVisControl->mVisSettings);
        }
      }

      snapshotsCompleted.erase(mSimulationState->GetTimeStepsPassed());
    }

    if (networkImagesCompleted.count(mSimulationState->GetTimeStepsPassed()) > 0)
    {
      for (std::multimap<unsigned long, unsigned long>::const_iterator it =
          networkImagesCompleted.find(mSimulationState->GetTimeStepsPassed()); it
          != networkImagesCompleted.end() && it->first == mSimulationState->GetTimeStepsPassed(); ++it)
      {
        if (hemelb::topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
        {

          const hemelb::vis::PixelSet<hemelb::vis::ResultPixel>* result =
              mVisControl->GetResult(it->second);

          if (steeringCpt->updatedMouseCoords)
          {
            float density, stress;

            if (mVisControl->MouseIsOverPixel(result, &density, &stress))
            {
              double mouse_pressure, mouse_stress;
              mLbm->CalculateMouseFlowField(density,
                                            stress,
                                            mouse_pressure,
                                            mouse_stress,
                                            mVisControl->mDomainStats.density_threshold_min,
                                            mVisControl->mDomainStats.density_threshold_minmax_inv,
                                            mVisControl->mDomainStats.stress_threshold_max_inv);

              mVisControl->SetMouseParams(mouse_pressure, mouse_stress);
            }
            steeringCpt->updatedMouseCoords = false;
          }

          imageSendCpt->DoWork(result);

        }
      }

      networkImagesCompleted.erase(mSimulationState->GetTimeStepsPassed());
    }

    double lPreSnapshotTime = hemelb::util::myClock();

    if (mSimulationState->GetTimeStep() % snapshots_period == 0)
    {
      char snapshot_filename[255];
      snprintf(snapshot_filename, 255, "snapshot_%06li.dat", mSimulationState->GetTimeStep());

      mSnapshotsWritten++;
      mLbm->WriteConfigParallel(stability, snapshotDirectory + std::string(snapshot_filename));
    }

    mSnapshotTime += (hemelb::util::myClock() - lPreSnapshotTime);

    if (stability == hemelb::lb::StableAndConverged)
    {
      is_finished = true;
      break;
    }
    if (mSimulationState->GetIsTerminating())
    {
      is_finished = true;
      break;
    }
    if (mSimulationState->GetTimeStepsPerCycle() > 400000)
    {
      is_unstable = true;
      break;
    }

    if (mSimulationState->GetTimeStep() == mSimulationState->GetTimeStepsPerCycle()
        && hemelb::topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
    {
      fprintf(mTimingsFile, "cycle id: %li\n", mSimulationState->GetCycleId());

      hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("cycle id: %li",
                                                                          mSimulationState->GetCycleId());

      fflush( NULL);
    }
  }

  PostSimulation(total_time_steps, hemelb::util::myClock() - simulation_time, is_unstable);

  hemelb::log::Logger::Log<hemelb::log::Warning, hemelb::log::Singleton>("Finish running simulation.");
}

/**
 * Called on error to abort the simulation and pull-down the MPI environment.
 */
void SimulationMaster::Abort()
{
  int err = MPI_Abort(MPI_COMM_WORLD, 1);

  // This gives us something to work from when we have an error - we get the rank
  // that calls abort, and we get a stack-trace from the exception having been thrown.
  hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::OnePerCore>("Aborting");
  exit(1);
}

/**
 * Steps that are taken when the simulation is complete.
 *
 * This function writes several bits of timing data to
 * the timing file.
 */
void SimulationMaster::PostSimulation(int iTotalTimeSteps, double iSimulationTime, bool iIsUnstable)
{
  if (hemelb::topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
  {
    fprintf(mTimingsFile, "\n");
    fprintf(mTimingsFile,
            "threads: %i, machines checked: %i\n\n",
            hemelb::topology::NetworkTopology::Instance()->GetProcessorCount(),
            hemelb::topology::NetworkTopology::Instance()->GetMachineCount());
    fprintf(mTimingsFile,
            "topology depths checked: %i\n\n",
            hemelb::topology::NetworkTopology::Instance()->GetDepths());
    fprintf(mTimingsFile, "fluid sites: %li\n\n", mLbm->TotalFluidSiteCount());
    fprintf(mTimingsFile,
            "cycles and total time steps: %li, %i \n\n",
             (mSimulationState->GetCycleId() - 1), // Note that the cycle-id is 1-indexed.
            iTotalTimeSteps);
    fprintf(mTimingsFile, "time steps per second: %.3f\n\n", iTotalTimeSteps / iSimulationTime);
  }

  if (iIsUnstable)
  {
    if (hemelb::topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
    {
      fprintf(mTimingsFile,
              "Attention: simulation unstable with %li timesteps/cycle\n",
              (unsigned long) mSimulationState->GetTimeStepsPerCycle());
      fprintf(mTimingsFile, "Simulation terminated\n");
    }
  }
  else
  {
    if (hemelb::topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
    {

      fprintf(mTimingsFile,
              "time steps per cycle: %li\n",
              (unsigned long) mSimulationState->GetTimeStepsPerCycle());
      fprintf(mTimingsFile, "\n");

      fprintf(mTimingsFile, "\n");

      fprintf(mTimingsFile, "\n");
      fprintf(mTimingsFile, "decomposition optimisation time (s):       %.3f\n", mDomainDecompTime);
      fprintf(mTimingsFile, "pre-processing buffer management time (s): %.3f\n", mNetInitialiseTime);
      fprintf(mTimingsFile, "input configuration reading time (s):      %.3f\n", mFileReadTime);

      fprintf(mTimingsFile,
              "total time (s):                            %.3f\n\n",
               (hemelb::util::myClock() - mCreationTime));

      fprintf(mTimingsFile, "Sub-domains info:\n\n");

      for (hemelb::proc_t n = 0; n
          < hemelb::topology::NetworkTopology::Instance()->GetProcessorCount(); n++)
      {
        fprintf(mTimingsFile,
                "rank: %lu, fluid sites: %lu\n",
                (unsigned long) n,
                (unsigned long) hemelb::topology::NetworkTopology::Instance()->FluidSitesOnEachProcessor[n]);
      }
    }
  }

  PrintTimingData();

}

/**
 * Outputs a breakdown of the simulation time spent on different activities.
 */
void SimulationMaster::PrintTimingData()
{
  // Note that CycleId is 1-indexed and will have just been incremented when we finish.
  double cycles = hemelb::util::NumericalFunctions::max(1.0,
                                                        (double) (mSimulationState->GetCycleId()
                                                            - 1));

  double lTimings[5] = { mLbm->GetTimeSpent() / cycles, mMPISendTime / cycles, mMPIWaitTime
      / cycles, mVisControl->GetTimeSpent()
      / hemelb::util::NumericalFunctions::max(1.0, (double) mImagesWritten), mSnapshotTime
      / hemelb::util::NumericalFunctions::max(1.0, (double) mSnapshotsWritten) };
  std::string lNames[5] = { "LBM", "MPISend", "MPIWait", "Images", "Snaps" };

  double lMins[5];
  double lMaxes[5];
  double lMeans[5];

  MPI_Reduce(lTimings, lMaxes, 5, hemelb::MpiDataType(lMaxes[0]), MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(lTimings, lMeans, 5, hemelb::MpiDataType(lMeans[0]), MPI_SUM, 0, MPI_COMM_WORLD);

  // Change the values for LBM and MPI on process 0 so they don't interfere with the min
  // operation (previously values were 0.0 so they won't affect max / mean
  // calc).

  MPI_Reduce(lTimings, lMins, 5, hemelb::MpiDataType(lMins[0]), MPI_MIN, 0, MPI_COMM_WORLD);

  if (hemelb::topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
  {
    for (int ii = 0; ii < 5; ii++)
      lMeans[ii] /= (double) (hemelb::topology::NetworkTopology::Instance()->GetProcessorCount());

    fprintf(mTimingsFile,
            "\n\nPer-proc timing data (secs per [cycle,cycle,cycle,image,snapshot]): \n\n");
    fprintf(mTimingsFile, "\t\tMin \tMean \tMax\n");
    for (int ii = 0; ii < 5; ii++)
    {
      fprintf(mTimingsFile,
              "%s\t\t%.3g\t%.3g\t%.3g\n",
              lNames[ii].c_str(),
              lMins[ii],
              lMeans[ii],
              lMaxes[ii]);
    }
  }
}
