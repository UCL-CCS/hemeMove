#include "lb/boundaries/BoundaryValues.h"
#include "util/utilityFunctions.h"
#include "util/fileutils.h"
#include <algorithm>
#include <fstream>

namespace hemelb
{
  namespace lb
  {
    namespace boundaries
    {

      BoundaryValues::BoundaryValues(geometry::SiteType IOtype,
                                     geometry::LatticeData* iLatDat,
                                     std::vector<iolets::InOutLet*> &iiolets,
                                     SimulationState* iSimState,
                                     util::UnitConverter* units) :
          net::IteratedAction(), mState(iSimState), mUnits(units)
      {
        nTotIOlets = (int) iiolets.size();

        std::vector<int> *procsList = new std::vector<int>[nTotIOlets];

        nIOlets = 0;

        // Determine which iolets need comms and create them
        for (int i = 0; i < nTotIOlets; i++)
        {
          // First create a copy of all iolets
          iolets::InOutLet* iolet = iiolets[i]->Clone();

          iolet->Initialise(mUnits);

          iolets.push_back(iolet);

          bool IOletOnThisProc = IsIOletOnThisProc(IOtype, iLatDat, i);
          procsList[i] = GatherProcList(IOletOnThisProc);

          // With information on whether a proc has an IOlet and the list of procs for each IOlte
          // on the BC task we can create the comms
          if (IOletOnThisProc || IsCurrentProcTheBCProc())
          {
            nIOlets++;
            ioletIDs.push_back(i);
            if (iolets[i]->DoComms())
            {
              mComms.push_back(new BoundaryComms(mState, procsList[i], IOletOnThisProc));
            }
            else
            {
              // NULL values fill up space to make indexing easier
              mComms.push_back(NULL);
            }
          }
        }

        densityCycle = new std::vector<distribn_t>[nIOlets];

        // Send out initial values
        Reset();

        // Clear up
        delete[] procsList;

      }

      BoundaryValues::~BoundaryValues()
      {

        delete[] densityCycle;

        for (int i = 0; i < nTotIOlets; i++)
        {
          delete iolets[i];
        }

        for (int i = 0; i < nIOlets; i++)
        {
          if (mComms[i] != NULL)
          {
            delete mComms[i];
          }
        }
      }

      bool BoundaryValues::IsIOletOnThisProc(geometry::SiteType IOtype, geometry::LatticeData* iLatDat, int iBoundaryId)
      {
        for (site_t i = 0; i < iLatDat->GetLocalFluidSiteCount(); i++)
        {
          const geometry::Site site = iLatDat->GetSite(i);

          if (site.GetSiteType() == IOtype && site.GetBoundaryId() == iBoundaryId)
          {
            return true;
          }
        }

        return false;
      }

      std::vector<int> BoundaryValues::GatherProcList(bool hasBoundary)
      {
        std::vector<int> procsList(0);

        // This is where the info about whether a proc contains the given inlet/outlet is sent
        // If it does contain the given inlet/outlet it sends a true value, else it sends a false.
        int IOletOnThisProc = hasBoundary; // true if inlet i is on this proc

        // These should be bool, but MPI only supports MPI_INT
        // For each inlet/outlet there is an array of length equal to total number of procs.
        // Each stores true/false value. True if proc of rank equal to the index contains
        // the given inlet/outlet.
        int nTotProcs = topology::NetworkTopology::Instance()->GetProcessorCount();
        int *boolList = new int[nTotProcs];

        MPI_Gather(&IOletOnThisProc,
                   1,
                   hemelb::MpiDataType(IOletOnThisProc),
                   boolList,
                   1,
                   hemelb::MpiDataType(boolList[0]),
                   GetBCProcRank(),
                   MPI_COMM_WORLD);

        if (IsCurrentProcTheBCProc())
        {
          // Now we have an array for each IOlet with true (1) at indices corresponding to
          // processes that are members of that group. We have to convert this into arrays
          // of ints which store a list of processor ranks.
          for (int j = 0; j < nTotProcs; j++)
          {
            if (boolList[j])
            {
              procsList.push_back(j);
            }
          }
        }

        // Clear up
        delete[] boolList;

        return procsList;
      }

      bool BoundaryValues::IsCurrentProcTheBCProc()
      {
        return topology::NetworkTopology::Instance()->GetLocalRank() == GetBCProcRank();
      }

      proc_t BoundaryValues::GetBCProcRank()
      {
        return 0;
      }

      void BoundaryValues::RequestComms()
      {
        for (int i = 0; i < nIOlets; i++)
        {
          unsigned long time_step = (mState->Get0IndexedTimeStep()) % densityCycle[i].size();
          SetDensityCycle(iolets[ioletIDs[i]],i,time_step);
        }
      }

      void BoundaryValues::SetDensityCycle(iolets::InOutLet* iolet, unsigned int iolet_index, unsigned long time_step)
      {
        if (IsCurrentProcTheBCProc())
        {

          if (iolet->DoComms())
          {
            iolet->UpdateCycle(densityCycle[iolet_index], mState);
            mComms[iolet_index]->Send(&densityCycle[iolet_index][time_step]);
          }
        }
        if (iolet->DoComms())
        {
          double incoming_density;
          mComms[iolet_index]->Receive(&incoming_density);
          iolet->SetDensity(incoming_density);
        }
        else
        {
          unsigned long time_step = (mState->Get0IndexedTimeStep()) % iolet->GetUpdatePeriod(mState->GetTotalTimeSteps());
          iolet->UpdateCycle(densityCycle[iolet_index], mState);
          iolet->SetDensity(densityCycle[iolet_index][time_step]);
        }
      }

      void BoundaryValues::EndIteration()
      {
        for (int i = 0; i < nIOlets; i++)
        {
          if (iolets[ioletIDs[i]]->DoComms())
          {
            mComms[i]->FinishSend();
          }
        }
      }

      void BoundaryValues::FinishReceive()
      {
        for (int i = 0; i < nIOlets; i++)
        {
          if (iolets[ioletIDs[i]]->DoComms())
          {
            mComms[i]->Wait();
          }
        }
      }

      void BoundaryValues::Reset()
      {
        for (int i = 0; i < nIOlets; i++)
        {
          iolets::InOutLet* iolet=iolets[ioletIDs[i]];
          densityCycle[i].resize(iolet->GetUpdatePeriod(mState->GetTotalTimeSteps()));
          SetDensityCycle(iolet,i,0);
        }

        for (int i = 0; i < nIOlets; i++)
        {
          if (iolets[ioletIDs[i]]->DoComms())
          {
            mComms[i]->WaitAllComms();
          }
        }
      }

// This assumes the program has already waited for comms to finish before
      distribn_t BoundaryValues::GetBoundaryDensity(const int index)
      {
        return iolets[index]->GetDensity();
      }

      distribn_t BoundaryValues::GetDensityMin(int iBoundaryId)
      {
        return iolets[iBoundaryId]->GetDensityMin();
      }

      distribn_t BoundaryValues::GetDensityMax(int iBoundaryId)
      {
        return iolets[iBoundaryId]->GetDensityMax();
      }

    }
  }
}
