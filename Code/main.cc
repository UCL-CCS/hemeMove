// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "net/mpi.h"
#include "net/IOCommunicator.h"
#include "configuration/CommandLine.h"
#include "configuration/SimBuilder.h"
#include "io/ensure_hexfloat.h"
#include "debug.h"

int main(int argc, char *argv[])
{
  // main function needed to perform the entire simulation. Some
  // simulation parameters and performance statistics are output on
  // standard output
  using namespace hemelb;
  // Bring up MPI
  net::MpiEnvironment mpi(argc, argv);
  log::Logger::Init();
  try
  {
    net::MpiCommunicator commWorld = net::MpiCommunicator::World();

    net::IOCommunicator hemelbCommunicator(commWorld);
    log::Logger::Log<log::Info, log::Singleton>(
      "Initialised MPI: %d nodes x %d processes per node = %d total",
      hemelbCommunicator.GetLeadersComm().Size(),
      hemelbCommunicator.GetNodeComm().Size(),
      hemelbCommunicator.Size()
    );
    io::GlobalHexFloatLocale ensure_hexfloat;

    try {
      // Parse command line
      auto options = configuration::CommandLine(argc, argv);

      // Start the debugger (if requested)
      debug::Init(options.GetDebug(), argv[0], commWorld);

      // Prepare main simulation object...
      auto controller = configuration::SimBuilder::CreateSim<Traits<>>(options, hemelbCommunicator);

      // ..and run it.
      controller->RunSimulation();
    }

    // Interpose this catch to print usage before propagating the error.
    catch (configuration::CommandLine::OptionError& e)
    {
      log::Logger::Log<log::Critical, log::Singleton>(configuration::CommandLine::GetUsage());
      throw;
    }
  }
  catch (std::exception& e)
  {
    log::Logger::Log<log::Critical, log::OnePerCore>(e.what());
    mpi.Abort(-1);
  }
  // MPI gets finalised by MpiEnv's d'tor.
}
