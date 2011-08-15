#include "lb/BoundaryComms.h"
#include "topology/NetworkTopology.h"
#include "util/utilityFunctions.h"
#include <math.h>

namespace hemelb
{
  namespace lb
  {

    BoundaryComms::BoundaryComms(SimulationState* iSimState, int iTotIOlets) :
      net::IteratedAction(), nTotIOlets(iTotIOlets), mState(iSimState)
    {
      // Work out stuff for simulation (ie. should give same result on all procs)
      proc_t BCrank = 0;

      if (topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
        BCrank = topology::NetworkTopology::Instance()->GetLocalRank();

      // Since only one proc will update BCrank, the sum of all BCrank is the BCproc
      MPI_Allreduce(&BCrank, &BCproc, 1, hemelb::MpiDataType(BCrank), MPI_SUM, MPI_COMM_WORLD);

    }

    BoundaryComms::~BoundaryComms()
    {

      if (topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
      {
        for (int i = 0; i < nTotIOlets; i++)
          delete[] procsList[i];

        delete[] requestOffset;
        delete[] procsList;
        delete[] nProcs;
      }

      // Communicators and groups
      delete[] request;
      delete[] status;
    }

    void BoundaryComms::Initialise(geometry::LatticeData::SiteType IOtype,
                                   geometry::LatticeData* iLatDat,
                                   std::vector<distribn_t>* iDensityCycleVector)
    {
      density_cycle_vector = iDensityCycleVector;
      density_cycle = new distribn_t[density_cycle_vector->size()];

      for (unsigned int i = 0; i < density_cycle_vector->size(); i++)
        density_cycle[i] = (*density_cycle_vector)[i];

      // Work out which and how many inlets/outlets on this process
      nIOlets = 0;
      IOlets = std::vector<int>(0);

      if (topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
      {
        nIOlets = nTotIOlets;

        nProcs = new int[nTotIOlets];
        procsList = new int*[nTotIOlets];

        for (int i = 0; i < nTotIOlets; i++)
          nProcs[i] = 0;
      }

      for (site_t i = 0; i < iLatDat->GetLocalFluidSiteCount(); i++)
      {
        if (iLatDat->GetSiteType(i) == IOtype
            && !util::VectorFunctions::member(IOlets, iLatDat->GetBoundaryId(i)))
        {
          nIOlets++;
          IOlets.push_back(iLatDat->GetBoundaryId(i));
        }
      }

      // Now BC process must find out which process belongs to what group

      if (topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
      {
        // These should be bool, but MPI only supports MPI_INT
        // For each inlet/outlet there is an array of length equal to total number of procs.
        // Each stores true/false value. True if proc of rank equal to the index contains
        // the given inlet/outlet.
        int **boolList = new int*[nTotIOlets];

        int nTotProcs = topology::NetworkTopology::Instance()->GetProcessorCount();

        for (int i = 0; i < nTotIOlets; i++)
        {
          boolList[i] = new int[nTotProcs];
        }

        MPI_Status tempStat;

        // Non-BC procs will be sending at this point
        for (int i = 0; i < nTotIOlets; i++)
        {
          for (int proc = 0; proc < nTotProcs; proc++)
          {
            if (proc != BCproc)
            {
              MPI_Recv(&boolList[i][proc], 1, MPI_INT, proc, 100, MPI_COMM_WORLD, &tempStat);
            }
            else
            {
              // All of them should be false, but in case we ever want BCproc to also have some sites on board
              boolList[i][proc] = util::VectorFunctions::member(IOlets, i);
            }
          }
        }

        // Now we have an array for each IOlet with true (1) at indices corresponding to
        // processes that are members of that group. We have to convert this into arrays
        // of ints which store a list of processor ranks.

        for (int i = 0; i < nTotIOlets; i++)
        {
          for (int j = 0; j < nTotProcs; j++)
          {
            if (boolList[i][j])
              nProcs[i]++;
          }

          procsList[i] = new int[nProcs[i]];

          int memberIndex = 0;

          for (int j = 0; j < nTotProcs; j++)
          {
            if (boolList[i][j])
            {
              procsList[i][memberIndex] = j;
              memberIndex++;
            }
          }
        }

        // Clear up
        for (int i = 0; i < nTotIOlets; i++)
        {
          delete[] boolList[i];
        }

        delete[] boolList;

      }
      else
      {
        // This is where the info about whether a proc contains a given inlet/outlet is sent
        // If it does contain the given inlet/outlet it sends a true value, else it sends a false.
        for (int i = 0; i < nTotIOlets; i++)
        {
          int IOletOnThisProc = util::VectorFunctions::member(IOlets, i); // true if inlet i is on this proc

          MPI_Ssend(&IOletOnThisProc, 1, MPI_INT, BCproc, 100, MPI_COMM_WORLD);
        }
      }

      if (topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
      {
        density = density_cycle;

        int nRequests = 0;

        requestOffset = new int[nTotIOlets];

        for (int i = 0; i < nTotIOlets; i++)
        {
          requestOffset[i] = nRequests;
          nRequests += nProcs[i];
        }

        request = new MPI_Request[nRequests];
        status = new MPI_Status[nRequests];
      }
      else
      {
        density = new distribn_t[nTotIOlets];
        request = new MPI_Request[nTotIOlets];
        status = new MPI_Status[nTotIOlets];
      }

      RequestComms();
      WaitAllComms();
    }

    // This assumes that the BCproc never receives
    distribn_t BoundaryComms::GetBoundaryDensity(const int index)
    {
      MPI_Wait(&request[index], &status[index]);

      return density[index];
    }

    void BoundaryComms::WaitAllComms()
    {
      if (topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
      {
        for (int i = 0; i < nTotIOlets; i++)
        {
          MPI_Waitall(nProcs[i], &request[requestOffset[i]], &status[requestOffset[i]]);
        }
      }
      else
      {
        for (int i = 0; i < nIOlets; i++)
        {
          MPI_Wait(&request[IOlets[i]], &status[IOlets[i]]);
        }
      }
    }

    void BoundaryComms::RequestComms()
    {
      if (topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
      {
        unsigned long time_step = mState->GetTimeStep() % mState->GetTimeStepsPerCycle();
        density = &density_cycle[time_step * nTotIOlets];

        int message = 0;

        for (int i = 0; i < nTotIOlets; i++)
        {
          for (int proc = 0; proc < nProcs[i]; proc++)
          {
            MPI_Isend(&density[i],
                      1,
                      hemelb::MpiDataType(density[0]),
                      procsList[i][proc],
                      100,
                      MPI_COMM_WORLD,
                      &request[message++]);
          }
        }
      }
      else
      {
        for (int i = 0; i < nIOlets; i++)
        {
          MPI_Irecv(&density[IOlets[i]],
                    1,
                    hemelb::MpiDataType(density[0]),
                    BCproc,
                    100,
                    MPI_COMM_WORLD,
                    &request[IOlets[i]]);
        }
      }
    }

    void BoundaryComms::EndIteration()
    {
      // Don't move on to next step with BC proc until all messages have been sent
      // Precautionary measure to make sure proc doesn't overwrite, before message is sent
      if (topology::NetworkTopology::Instance()->IsCurrentProcTheIOProc())
      {
        for (int i = 0; i < nTotIOlets; i++)
          MPI_Waitall(nProcs[i], &request[requestOffset[i]], &status[requestOffset[i]]);
      }
    }

    void BoundaryComms::Reset()
    {
      density_cycle = new distribn_t[density_cycle_vector->size()];

      for (unsigned int i = 0; i < density_cycle_vector->size(); i++)
        density_cycle[i] = (*density_cycle_vector)[i];

      // density_cycle should be reset by now
      RequestComms();
      WaitAllComms();
    }

  }
}
