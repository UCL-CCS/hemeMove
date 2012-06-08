#include "configuration/CommandLine.h"
#include "SimulationMaster.h"

int main(int argc, char *argv[])
{
  // main function needed to perform the entire simulation. Some
  // simulation paramenters and performance statistics are outputted on
  // standard output

  hemelb::configuration::CommandLine options = hemelb::configuration::CommandLine(argc,argv);
  SimulationMaster lMaster = SimulationMaster(options);

  lMaster.RunSimulation();

  return (0);
}

/*#include "configuration/CommandLine.h"
//#include "SimulationMaster.h"
#include <multiscale/MultiscaleSimulationMaster.h>
#include "unittests/multiscale/MockMPWide.h"
#include "unittests/multiscale/MockIntercommunicand.h"
#include <resources/Resource.h>
#include "unittests/multiscale/MPWideIntercommunicator.h"

int main(int argc, char *argv[])
{
  // main function needed to perform the entire simulation. Some
  // simulation paramenters and performance statistics are outputted on
  // standard output

  hemelb::configuration::CommandLine options = hemelb::configuration::CommandLine(argc,argv);
  
  MultiscaleSimulationMaster<MPWideIntercommunicator> *lMaster;

  pbuffer = new std::map<std::string, double>();
  std::map<std::string, double> &buffer = *pbuffer;

  orchestrationLB=new std::map<std::string,bool>();
  std::map<std::string,bool> &rorchestrationLB=*orchestrationLB;
  rorchestrationLB["boundary1_pressure"] = false;
  rorchestrationLB["boundary2_pressure"] = false;
  rorchestrationLB["boundary1_velocity"] = true;
  rorchestrationLB["boundary2_velocity"] = true;

  MPWideIntercommunicator intercomms(*pbuffer,*orchestrationLB);

  //TODO: Add an IntercommunicatorImplementation?

  lMaster = new MultiscaleSimulationMaster<MPWideIntercommunicator>(options, intercomms);

  lMaster.RunSimulation();

  return (0);
}*/
