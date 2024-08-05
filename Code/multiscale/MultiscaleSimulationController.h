// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_MULTISCALE_MULTISCALESIMULATIONCONTROLLER_H
#define HEMELB_MULTISCALE_MULTISCALESIMULATIONCONTROLLER_H

#include <vector>
#include "lb/iolets/InOutLetMultiscale.h"
#include "multiscale/Intercommunicator.h"
#include "SimulationController.h"
#include "util/span.h"

namespace hemelb::multiscale
{
    /***
     * Instead of adding multiscale functionality to the standard simulation controller, we keep this here,
     * so the main code can be read without thinking about multiscale.
     */
    template<class Intercommunicator>
    class MultiscaleSimulationController : public SimulationController
    {
    public:
        MultiscaleSimulationController(configuration::CommandLine &options,
                                   const net::IOCommunicator& ioComm,
                                   Intercommunicator & aintercomms) :
            SimulationController(options, ioComm), intercomms(aintercomms),
                multiscaleIoletType("inoutlet")
        {
          // We only have one shared object type so far, an iolet.
          lb::InOutLetMultiscale::DefineType(multiscaleIoletType);

          log::Logger::Log<log::Info, log::OnePerCore>("CONSTRUCTOR: inlet and outlet count: %d and %d",
                                                                               inletValues->GetLocalIoletCount(),
                                                                               outletValues->GetLocalIoletCount());
          log::Logger::Log<log::Warning, log::OnePerCore>("inlets: %d",
                                                                                  inletValues->GetLocalIolet(0)->IsCommsRequired(),
                                                                                  inletValues->GetLocalIolet(0)->GetDensityMax(),
                                                                                  inletValues->GetLocalIolet(0)->GetPressureMax());
          log::Logger::Log<log::Warning, log::OnePerCore>("outlets: %d",
                                                                                  outletValues->GetLocalIolet(0)->IsCommsRequired(),
                                                                                  outletValues->GetLocalIolet(0)->GetDensityMax(),
                                                                                  outletValues->GetLocalIolet(0)->GetPressureMax());

          // we only want to register those iolets which are needed on this process.
          // Fortunately, the BoundaryValues instance has worked this out for us.
          for (unsigned int i = 0; i < inletValues->GetLocalIoletCount(); i++)
          {
            // could be a if dynamic_cast<> rather than using a castable? virtual method pattern, if we prefer.
            if (inletValues->GetLocalIolet(i)->IsRegistrationRequired())
            {
              static_cast<lb::InOutLetMultiscale*>(inletValues->GetLocalIolet(i))->Register(intercomms,
                                                                                                    multiscaleIoletType);
            }
          }
          for (unsigned int i = 0; i < outletValues->GetLocalIoletCount(); i++)
          {
            if (outletValues->GetLocalIolet(i)->IsRegistrationRequired())
            {
              static_cast<lb::InOutLetMultiscale*>(outletValues->GetLocalIolet(i))->Register(intercomms,
                                                                                                     multiscaleIoletType);
            }
          }
          // We only have one shared object type so far, an iolet.
          lb::InOutLetMultiscale::DefineType(multiscaleIoletType);

          /* Process 0 has a list of all the Iolets. The count of all this is highly useful to pre-size all the
           * needed arrays later on, so we are broadcasting this to all the other processes. */
          std::vector<unsigned> GlobalIoletCount;
          GlobalIoletCount.push_back(inletValues->GetLocalIoletCount());
          GlobalIoletCount.push_back(outletValues->GetLocalIoletCount());
          ioComms.Broadcast(to_span(GlobalIoletCount), 0);

          std::vector<std::vector<site_t> > invertedInletBoundaryList(GlobalIoletCount[0]);
          std::vector<std::vector<site_t> > invertedOutletBoundaryList(GlobalIoletCount[1]);

          log::Logger::Log<log::Debug, log::OnePerCore>("inlets start %i/%i",
                                                                                inletValues->GetLocalIoletCount(),
                                                                                GlobalIoletCount[0]);
          log::Logger::Log<log::Debug, log::OnePerCore>("outlets start %i/%i",
                                                                                outletValues->GetLocalIoletCount(),
                                                                                GlobalIoletCount[1]);

          //TODO: Throw a warning when process 0 count mismatches with the aggregate of the others.
          bool velocity = false;
          if (velocity == true)
          {
            /* Do not include the non-iolet adjacent sites (resp. MidFluid and Wall-adjacent). */

            /* Do include iolet adjacent sites (inlet) */
            invertedInletBoundaryList = PopulateInvertedBoundaryList(domainData.get(),
                                                                     invertedInletBoundaryList,
                                                                     domainData->GetMidDomainSiteRange(2));

            /* Do include iolet adjacent sites (outlet) */
            invertedOutletBoundaryList = PopulateInvertedBoundaryList(domainData.get(),
                                                                      invertedOutletBoundaryList,
                                                                      domainData->GetMidDomainSiteRange(3));

            /* Do include iolet adjacent sites (inlet-wall) */
            invertedInletBoundaryList = PopulateInvertedBoundaryList(domainData.get(),
                                                                     invertedInletBoundaryList,
                                                                     domainData->GetMidDomainSiteRange(4));

            /* Do include iolet adjacent sites (outlet-wall) */
            invertedOutletBoundaryList = PopulateInvertedBoundaryList(domainData.get(),
                                                                      invertedOutletBoundaryList,
                                                                      domainData->GetMidDomainSiteRange(5));

            log::Logger::Log<log::Debug, log::OnePerCore>("Populated inlets (numinlets/sizeinlet0): %i/%i",
                                                                                  invertedInletBoundaryList.size(),
                                                                                  invertedInletBoundaryList[0].size());
            log::Logger::Log<log::Debug, log::OnePerCore>("Populated outlets (numoutlets/sizeoutlet0): %i/%i",
                                                                                  invertedOutletBoundaryList.size(),
                                                                                  invertedOutletBoundaryList[0].size());
          }

          //TODO: Debug
          //invertedInletBoundaryList =
          //    ExchangeAndCompleteInverseBoundaryList(invertedInletBoundaryList);
          //invertedOutletBoundaryList =
          //    ExchangeAndCompleteInverseBoundaryList(invertedOutletBoundaryList);

          //lb::MacroscopicPropertyCache& propertyCache =
          //    latticeBoltzmannModel->GetPropertyCache();

          //WORKAROUND: Reserve space in the velocityCache to ensure that all elements are properly allocated.
          //If this is not set, then velocityCache will contain 0 elements at the start even though the
          //MacroscopicPropertyCache reports a number of site count higher than 0.
          //propertyCache->velocityCache.Reserve(propertyCache->GetSiteCount());

          // we only want to register those iolets which are needed on this process.
          // Fortunately, the BoundaryValues instance has worked this out for us.
          for (unsigned int i = 0; i < inletValues->GetLocalIoletCount(); i++)
          {
            log::Logger::Log<log::Debug, log::OnePerCore>("1) %i %i",
                                                                                  i,
                                                                                  GlobalIoletCount[0]);
            // could be a if dynamic_cast<> rather than using a castable? virtual method pattern, if we prefer.
            if (inletValues->GetLocalIolet(i)->IsRegistrationRequired())
            {
              log::Logger::Log<log::Debug, log::OnePerCore>("2) inlets: %i %i %i",
                                                                                    invertedInletBoundaryList.size(),
                                                                                    invertedInletBoundaryList[0].size(),
                                                                                    i);
              static_cast<lb::InOutLetMultiscale*>(inletValues->GetLocalIolet(i))->Register(intercomms,
                                                                                                    multiscaleIoletType);
              log::Logger::Log<log::Debug, log::OnePerCore>("3) inlets: %i %i",
                                                                                    invertedInletBoundaryList.size(),
                                                                                    invertedInletBoundaryList[i].size());
              /*static_cast<lb::InOutLetVelocityAware*>(inletValues->GetLocalIolet(i))->InitialiseNeighbouringSites(neighbouringDataManager,
               m_fieldData,
               &propertyCache,
               invertedInletBoundaryList[i]);*/
            }
          }

          for (unsigned int i = 0; i < outletValues->GetLocalIoletCount(); i++)
          {
            log::Logger::Log<log::Debug, log::OnePerCore>("1) %i %i",
                                                                                  i,
                                                                                  GlobalIoletCount[1]);
            if (outletValues->GetLocalIolet(i)->IsRegistrationRequired())
            {
              log::Logger::Log<log::Debug, log::OnePerCore>("2) outlets: %i %i %i",
                                                                                    invertedOutletBoundaryList.size(),
                                                                                    invertedOutletBoundaryList[0].size(),
                                                                                    i);
              static_cast<lb::InOutLetMultiscale*>(outletValues->GetLocalIolet(i))->Register(intercomms,
                                                                                                     multiscaleIoletType);
              log::Logger::Log<log::Debug, log::OnePerCore>("3) outlets: %i %i",
                                                                                    invertedOutletBoundaryList.size(),
                                                                                    invertedOutletBoundaryList[i].size());
              /*static_cast<lb::InOutLetVelocityAware*>(outletValues->GetLocalIolet(i))->InitialiseNeighbouringSites(neighbouringDataManager,
               m_fieldData,
               &propertyCache,
               invertedOutletBoundaryList[i]);*/
            }
          }

          log::Logger::Log<log::Debug, log::OnePerCore>("MSController ShareICs started...");
          intercomms.ShareInitialConditions();
          log::Logger::Log<log::Debug, log::OnePerCore>("MSController Init finished!");
        }

        void PrintVectorList(std::vector<std::vector<site_t> > v)
        {
          log::Logger::Log<log::Info, log::OnePerCore>("Printing Vector List:");
          for (unsigned int i = 0; i < v.size(); i++)
          {
            for (unsigned int j = 0; j < v[i].size(); j++)
            {
              log::Logger::Log<log::Info, log::OnePerCore>("Boundary: %i % i %i",
                                                                                   i,
                                                                                   j,
                                                                                   v[i][j]);
            }
          }
        }

        void DoTimeStep() override
        {
          bool advance = intercomms.DoMultiscale(GetState().GetTime());
          log::Logger::Log<log::Info, log::Singleton>("At time step %i, should advance %i, time %f",
                                                                              GetState().GetTimeStep(),
                                                                              static_cast<int>(advance),
                                                                              GetState().GetTime());

          if (advance)
          {
            /* NOTE: Following triggers DoComms in each InOutLetMultiscale.
             * IoLetMS is aggressive with this, and DoComms will actually
             * complete all the communications, not just initiate them.
             * This is to prevent any inconsistent state in the coupling
             * (it's hard enough to get the physics right with a consistent
             * state ;)). */

            log::Logger::Log<log::Debug, log::OnePerCore>("inlet and outlet count: %d and %d",
                                                                                  inletValues->GetLocalIoletCount(),
                                                                                  outletValues->GetLocalIoletCount());
            log::Logger::Log<log::Debug, log::OnePerCore>("inlets: %d",
                                                                                  inletValues->GetLocalIolet(0)->IsCommsRequired(),
                                                                                  inletValues->GetLocalIolet(0)->GetDensityMax(),
                                                                                  inletValues->GetLocalIolet(0)->GetPressureMax());
            log::Logger::Log<log::Debug, log::OnePerCore>("outlets: %d",
                                                                                  outletValues->GetLocalIolet(0)->IsCommsRequired(),
                                                                                  outletValues->GetLocalIolet(0)->GetDensityMax(),
                                                                                  outletValues->GetLocalIolet(0)->GetPressureMax());

            SetCommsRequired(inletValues.get(), true);
            SetCommsRequired(outletValues.get(), true);

            inletValues->RequestComms();
            outletValues->RequestComms();
            SetCommsRequired(inletValues.get(), false);
            SetCommsRequired(outletValues.get(), false);

            for (unsigned int i = 0; i < inletValues->GetLocalIoletCount(); i++)
            {
              log::Logger::Log<log::Debug, log::OnePerCore>("Inlet[%i]: Measured Density is %f. Pressure is %f.",
                                                                                    i,
                                                                                    inletValues->GetLocalIolet(i)->GetDensity(GetState().GetTimeStep()),
                                                                                    inletValues->GetLocalIolet(i)->GetPressureMax());
            }
            for (unsigned int i = 0; i < outletValues->GetLocalIoletCount(); i++)
            {
              log::Logger::Log<log::Debug, log::OnePerCore>("Outlet[%i]: Measured Density is %f. Pressure is %f.",
                                                                                    i,
                                                                                    outletValues->GetLocalIolet(i)->GetDensity(GetState().GetTimeStep()),
                                                                                    outletValues->GetLocalIolet(i)->GetPressureMax());
            }

            /* Temporary Orchestration hardcode for testing 1/100 step ratio
             * TODO: Make an orchestration system for the multiscale coupling. */
            //for (int i = 0; i < 100; i++)
            //{
            log::Logger::Log<log::Debug, log::Singleton>("Step: HemeLB advanced to time %f.",
                                                                                 GetState().GetTime());
            SimulationController::DoTimeStep();
            //}
          }
          else
          {
            log::Logger::Log<log::Debug, log::Singleton>("HemeLB waiting pending multiscale siblings.");
            return;
          };
        }
        Intercommunicator &intercomms;
        typename Intercommunicator::IntercommunicandTypeT multiscaleIoletType;

    private:

        /* Loops over iolets to set the need for communications. */
        void SetCommsRequired(lb::BoundaryValues* ioletValues, bool b)
        {
          log::Logger::Log<log::Debug, log::OnePerCore>("Starting SetCommsRequired.");
          for (unsigned int i = 0; i < ioletValues->GetLocalIoletCount(); i++)
          {
            log::Logger::Log<log::Debug, log::OnePerCore>("In loop: %d",
                                                                                  ioletValues->GetLocalIoletCount());
            log::Logger::Log<log::Debug, log::OnePerCore>("A: iolet %d %d",
                                                                                  i,
                                                                                  (ioletValues->GetLocalIolet(i))->IsCommsRequired());
            log::Logger::Log<log::Debug, log::OnePerCore>("B: iolet %d %d",
                                                                                  i,
                                                                                  static_cast<lb::InOutLetMultiscale*>(ioletValues->GetLocalIolet(i))->IsCommsRequired());
            dynamic_cast<lb::InOutLetMultiscale*>(ioletValues->GetLocalIolet(i))->SetCommsRequired(b);
            log::Logger::Log<log::Debug, log::OnePerCore>("done with SetCommsRequired iteration.");

          }
          log::Logger::Log<log::Debug, log::OnePerCore>("Finishing SetCommsRequired.");
        }

        //Populate an invertedBoundaryList
        //out: iBL
        //in: domain_type, [Site Object]->SiteData->GetBoundaryID,
        //MPI_Gatherv if data is only this process.
        std::vector<std::vector<site_t> >
        PopulateInvertedBoundaryList(
            geometry::Domain* latticeData,
            std::vector<std::vector<site_t> > invertedBoundaryList, std::pair<site_t, site_t> i_range
        ) {
            for (auto i = i_range.first; i < i_range.second; ++i) {
                /* 1. Obtain Boundary ID number. */
                auto&& s = latticeData->GetSite(i).GetSiteData();
                int boundaryID = s.GetIoletId();

                /* 2. Grow the list to an appropriate size if needed. */
                while (std::ssize(invertedBoundaryList) <= boundaryID) {
              log::Logger::Log<log::Warning, log::OnePerCore>("WARNING: Growing the invertedBoundaryList, because we created in wrongly in MultiscaleSimulation Controller.");
              std::vector<site_t> a(0);
              invertedBoundaryList.push_back(a);
              if (invertedBoundaryList.size() > 100000)
              {
                log::Logger::Log<log::Warning, log::OnePerCore>("ERROR: invertedBoundaryList is growing to ridiculous proportions due to a faulty boundaryID.");
                exit(-1);
              }
                }

                /* 3. Insert this site in the inverted Boundary List. */
                invertedBoundaryList[boundaryID].push_back(
                        latticeData->GetGlobalNoncontiguousSiteIdFromGlobalCoords(
                                latticeData->GetSite(i).GetGlobalSiteCoords()
                        )
                );
          }
          return invertedBoundaryList;
        }

        std::vector<std::vector<site_t> > ExchangeAndCompleteInverseBoundaryList(
            std::vector<std::vector<site_t> > inList)
        {
          std::vector<std::vector<site_t> > outList;
          int *recvSizes = new int[ioComms.Size()];
          int *recvDispls = new int[ioComms.Size()];

          /* TODO: ASSUMPTION:
           * inList.size() is equal everywhere. This is not necessarily the case.
           * Use an AllReduce MAX and resize inList accordingly to make the remnant
           * of the code work here for unequal inList sizes. */

          for (unsigned int i = 0; i < inList.size(); i++)
          {
            int sendSize = ((int) inList[i].size());
            site_t *sendList = new site_t[inList[i].size()];
            for (unsigned int j = 0; j < inList[i].size(); j++)
            {
              sendList[j] = inList[i][j];
            }
            HEMELB_MPI_CALL(MPI_Allgather,
                            ( &sendSize, 1, MPI_INT, recvSizes, 1, MPI_INT, ioComms ));

            int64_t totalSize = 0;

            int np = ioComms.Size();
            int64_t offset = 0;

            for (int j = 0; j < np; j++)
            {
              totalSize += recvSizes[j];
              recvDispls[j] = offset;
              offset += recvSizes[j];
            }

            site_t *recvList = new site_t[totalSize]; //inList[i].size()

            HEMELB_MPI_CALL(MPI_Allgatherv,
                            ( sendList, inList[i].size(), MPI_LONG_LONG, recvList, recvSizes, recvDispls, MPI_LONG_LONG, ioComms ));

            std::vector<site_t> subList;
            for (int j = 0; j < totalSize; j++)
            {
              subList.push_back(recvList[j]);
            }
            outList.push_back(subList);
          }

          return outList;
        }
    };

}

#endif // HEMELB_MULTISCALE_MULTISCALESIMULATIONCONTROLLER_H
