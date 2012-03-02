#include "topology/NetworkTopology.h"
#include "mpiInclude.h"

namespace hemelb
{
  namespace topology
  {
    NetworkTopology NetworkTopology::instance;
    bool NetworkTopology::initialised = false;

    // Must be specified to prevent the default constructor being public.
    NetworkTopology::NetworkTopology()
    {

    }

    NetworkTopology* NetworkTopology::Instance()
    {
      return &instance;
    }

    void NetworkTopology::Init(int argCount, char ** argList, bool* oSuccess)
    {
      if (!initialised)
      {
        initialised = true;
        MPI_Init(&argCount, &argList);

        int tempSize = 0, tempRank = 0;
        MPI_Comm_size(MPI_COMM_WORLD, &tempSize);
        MPI_Comm_rank(MPI_COMM_WORLD, &tempRank);

        processorCount = (proc_t) tempSize;
        localRank = (proc_t) tempRank;

        *oSuccess = InitialiseMachineInfo();
      }
    }

    NetworkTopology::~NetworkTopology()
    {
      if (initialised)
      {
        MPI_Finalize();

        delete[] ProcCountOnEachMachine;
        delete[] MachineIdOfEachProc;
      }
    }

    bool NetworkTopology::IsCurrentProcTheIOProc() const
    {
      return localRank == 0;
    }

    proc_t NetworkTopology::GetLocalRank() const
    {
      return localRank;
    }

    proc_t NetworkTopology::GetProcessorCount() const
    {
      return processorCount;
    }

    int NetworkTopology::GetDepths() const
    {
      return depths;
    }

    unsigned int NetworkTopology::GetMachineCount() const
    {
      return machineCount;
    }

  }
}
