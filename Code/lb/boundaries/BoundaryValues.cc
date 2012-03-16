#include "lb/boundaries/BoundaryValues.h"
#include "lb/boundaries/BoundaryComms.h"
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
                                     const std::vector<iolets::InOutLet*> &iiolets,
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
            if (iolet->GetIsCommsRequired())
            {
              iolet->SetComms(new BoundaryComms(mState, procsList[i], IOletOnThisProc));
            }
          }
        }

        // Send out initial values
        Reset();

        // Clear up
        delete[] procsList;

      }

      BoundaryValues::~BoundaryValues()
      {

        for (int i = 0; i < nTotIOlets; i++)
        {
          delete iolets[i];
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
          HandleComms(iolets[ioletIDs[i]]);
        }
      }

      void BoundaryValues::HandleComms(iolets::InOutLet* iolet)
      {

        if (iolet->GetIsCommsRequired())
        {
          iolet->DoComms(IsCurrentProcTheBCProc());
        }

      }

      void BoundaryValues::EndIteration()
      {
        for (int i = 0; i < nIOlets; i++)
        {
          if (iolets[ioletIDs[i]]->GetIsCommsRequired())
          {
            iolets[ioletIDs[i]]->GetComms()->FinishSend();
          }
        }
      }

      void BoundaryValues::FinishReceive()
      {
        for (int i = 0; i < nIOlets; i++)
        {
          if (iolets[ioletIDs[i]]->GetIsCommsRequired())
          {
            iolets[ioletIDs[i]]->GetComms()->Wait();
          }
        }
      }

      void BoundaryValues::Reset()
      {
        for (int i = 0; i < nIOlets; i++)
        {
          iolets::InOutLet* iolet = iolets[ioletIDs[i]];
          iolet->Reset(*mState);
          if (iolet->GetIsCommsRequired())
          {
            iolet->GetComms()->WaitAllComms();

          }
        }
      }

      // This assumes the program has already waited for comms to finish before
      distribn_t BoundaryValues::GetBoundaryDensity(const int index)
      {
        return iolets[index]->GetDensity(mState->Get0IndexedTimeStep());
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
