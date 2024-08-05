// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "lb/iolets/BoundaryValues.h"

#include <algorithm>

#include "geometry/Domain.h"
#include "lb/iolets/BoundaryComms.h"
#include "log/Logger.h"
#include "util/numerical.h"

namespace hemelb::lb
{
      BoundaryValues::BoundaryValues(geometry::SiteType ioletType,
                                     geometry::Domain const& latticeData,
                                     const std::vector<IoletPtr> &incoming_iolets,
                                     SimulationState* simulationState,
                                     const net::MpiCommunicator& comms,
                                     const util::UnitConverter& unitConverter) :
              net::IteratedAction(), ioletType(ioletType),
              state(simulationState), bcComms(comms)
      {
          const auto totalIoletCount = incoming_iolets.size();
          std::vector<std::vector<int>> procsList(totalIoletCount);

        // Determine which iolets need comms and create them
        for (unsigned ioletIndex = 0; ioletIndex < totalIoletCount; ioletIndex++)
        {
          // First create a copy of all iolets
          auto iolet = incoming_iolets[ioletIndex].clone();

          iolet->Initialise(&unitConverter);

          bool isIoletOnThisProc = IsIoletOnThisProc(latticeData, ioletIndex);
          log::Logger::Log<log::Debug, log::OnePerCore>("BOUNDARYVALUES.CC - isioletonthisproc? : %d",
                                                                                isIoletOnThisProc);
          procsList[ioletIndex] = GatherProcList(isIoletOnThisProc);

          // With information on whether a proc has an iolet and the list of procs for each iolet
          // on the BC task we can create the comms
          if (isIoletOnThisProc || bcComms.IsCurrentProcTheBCProc())
          {
            localIoletIDs.push_back(ioletIndex);
//            hemelb::log::Logger::Log<hemelb::log::Warning, hemelb::log::OnePerCore>("BOUNDARYVALUES.H - ioletIndex: %d", ioletIndex);

//            if (iolet->IsCommsRequired()) //DEREK: POTENTIAL MULTISCALE ISSUE (this if-statement)
//            {
            iolet->SetComms(new BoundaryComms(state,
                                              procsList[ioletIndex],
                                              bcComms,
                                              isIoletOnThisProc));
//            }
          }
	  iolets.push_back(std::move(iolet));
        }

        // Send out initial values
        Reset();
      }

      bool BoundaryValues::IsIoletOnThisProc(geometry::Domain const& latticeData, int boundaryId)
      {
        for (site_t i = 0; i < latticeData.GetLocalFluidSiteCount(); i++)
        {
          auto&& site = latticeData.GetSite(i);

          if (site.GetSiteType() == ioletType && site.GetIoletId() == boundaryId)
          {
            return true;
          }
        }

        return false;
      }

      std::vector<int> BoundaryValues::GatherProcList(bool hasBoundary)
      {
        std::vector<int> processorsNeedingIoletList(0);

        // This is where the info about whether a proc contains the given inlet/outlet is sent
        // If it does contain the given inlet/outlet it sends a true value, else it sends a false.
        int isIoletOnThisProc = hasBoundary; // true if inlet i is on this proc

        // These should be bool, but MPI only supports MPI_INT
        // For each inlet/outlet there is an array of length equal to total number of procs.
        // Each stores true/false value. True if proc of rank equal to the index contains
        // the given inlet/outlet.

        std::vector<int> processorsNeedingIoletFlags = bcComms.Gather(isIoletOnThisProc,
                                                                      bcComms.GetBCProcRank());

        if (bcComms.IsCurrentProcTheBCProc())
        {
          // Now we have an array for each iolet with true (1) at indices corresponding to
          // processes that are members of that group. We have to convert this into arrays
          // of ints which store a list of processor ranks.
          for (proc_t process = 0; process < proc_t(processorsNeedingIoletFlags.size()); ++process)
          {
            if (processorsNeedingIoletFlags[process])
            {
              processorsNeedingIoletList.push_back(process);
            }
          }
        }

        return processorsNeedingIoletList; // return by copy
      }

      void BoundaryValues::RequestComms()
      {
        for (int i = 0; i < ssize(localIoletIDs); i++)
        {
          HandleComms(GetLocalIolet(i));
        }
      }

      void BoundaryValues::HandleComms(InOutLet* iolet)
      {

        if (iolet->IsCommsRequired())
        {
          iolet->DoComms(bcComms, state->GetTimeStep());
        }

      }

      void BoundaryValues::EndIteration()
      {
        for (int i = 0; i < ssize(localIoletIDs); i++)
        {
          if (GetLocalIolet(i)->IsCommsRequired())
          {
            GetLocalIolet(i)->GetComms()->FinishSend();
          }
        }
      }

      void BoundaryValues::FinishReceive()
      {
        for (int i = 0; i < ssize(localIoletIDs); i++)
        {
          if (GetLocalIolet(i)->IsCommsRequired())
          {
            GetLocalIolet(i)->GetComms()->Wait();
          }
        }
      }

      void BoundaryValues::Reset()
      {
        for (int i = 0; i < ssize(localIoletIDs); i++)
        {
          GetLocalIolet(i)->Reset(*state);
          if (GetLocalIolet(i)->IsCommsRequired())
          {
            GetLocalIolet(i)->GetComms()->WaitAllComms();

          }
        }
      }

      // This assumes the program has already waited for comms to finish before
      LatticeDensity BoundaryValues::GetBoundaryDensity(const int index)
      {
        return iolets[index]->GetDensity(state->GetTimeStep());
      }

      LatticeDensity BoundaryValues::GetDensityMin(int iBoundaryId)
      {
        return iolets[iBoundaryId]->GetDensityMin();
      }

      LatticeDensity BoundaryValues::GetDensityMax(int iBoundaryId)
      {
        return iolets[iBoundaryId]->GetDensityMax();
      }

}
