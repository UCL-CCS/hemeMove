#ifndef HEMELB_NET_H
#define HEMELB_NET_H

#include <vector>
#include <map>

#include "constants.h"
#include "mpiInclude.h"
#include "D3Q15.h"
#include "SimConfig.h"

#include "lb/GlobalLatticeData.h"
#include "lb/LocalLatticeData.h"
#include "topology/NetworkTopology.h"

namespace hemelb
{
  namespace net
  {

    class Net
    {
      public:
        Net(hemelb::topology::NetworkTopology * iTopology);
        ~Net();

        int* Initialise(hemelb::lb::GlobalLatticeData &iGlobLatDat,
                        hemelb::lb::LocalLatticeData* &bLocalLatDat);

        void Receive();
        void Send();
        void Wait(hemelb::lb::LocalLatticeData *bLocalLatDat);

        /**
         * Request that iCount entries of type T be included in the send to iToRank,
         * starting at oPointer.
         *
         * @param oPointer Pointer to first element to be sent.
         * @param iCount Number of element to be sent.
         * @param iToRank Rank to send to.
         */
        template<class T>
        void RequestSend(T* oPointer, int iCount, int iToRank)
        {
          if (sendReceivePrepped)
          {
            std::cerr
                << "Error: tried to add send-data after the datatype was already constructed. This is a bug.\n";
            exit(0);
          }

          ProcComms *lComms = GetProcComms(iToRank);

          AddToList(oPointer, iCount, lComms->SendData);
        }

        /**
         * Request that iCount entries of type T be included in the receive from iFromRank
         * starting at oPointer.
         *
         * @param oPointer Pointer to the first element of the array to receive into.
         * @param iCount Number of elements to receive into the array.
         * @param iFromRank Rank to receive from.
         */
        template<class T>
        void RequestReceive(T* oPointer, int iCount, int iFromRank)
        {
          if (sendReceivePrepped)
          {
            std::cerr
                << "Error: tried to add receive-data after the datatype was already constructed. This is a bug.\n";
            exit(0);
          }

          ProcComms *lComms = GetProcComms(iFromRank);

          AddToList(oPointer, iCount, lComms->ReceiveData);
        }

      private:
        /**
         * Struct representing all that's needed to successfully communicate with another processor.
         */
        struct ProcComms
        {
          public:
            struct MetaData
            {
                std::vector<void*> PointerList;
                std::vector<int> LengthList;
                std::vector<MPI_Datatype> TypeList;

                void clear()
                {
                  PointerList.clear();
                  LengthList.clear();
                  TypeList.clear();
                }
            };

            MetaData SendData;
            MetaData ReceiveData;

            MPI_Datatype SendType;
            MPI_Datatype ReceiveType;
        };

        ProcComms* GetProcComms(int iRank);

        void AddToList(int* iNew, int iLength, ProcComms::MetaData &bMetaData);

        void AddToList(double* iNew, int iLength, ProcComms::MetaData &bMetaData);

        void GetThisRankSiteData(const hemelb::lb::GlobalLatticeData & iGlobLatDat,
                                 unsigned int *& bThisRankSiteData);
        void InitialiseNeighbourLookup(hemelb::lb::LocalLatticeData *bLocalLatDat,
                                       short int **bSharedFLocationForEachProc,
                                       const unsigned int *iSiteDataForThisRank,
                                       const hemelb::lb::GlobalLatticeData & iGlobLatDat);
        void CountCollisionTypes(hemelb::lb::LocalLatticeData * bLocalLatDat,
                                 const hemelb::lb::GlobalLatticeData & iGlobLatDat,
                                 const unsigned int * lThisRankSiteData);

        void InitialisePointToPointComms(short int **& lSharedFLocationForEachProc);

        void EnsurePreparedToSendReceive();

        void CreateMPIType(const ProcComms::MetaData &iMetaData, MPI_Datatype &oNewDatatype);

        bool sendReceivePrepped;
        std::map<int, ProcComms*> mProcessorComms;
        hemelb::topology::NetworkTopology * mNetworkTopology;

        // Requests equal to twice the total number of processors in the communicator
        // are available for general communication within the Net object (both
        // initialisation and during each iteration).
        MPI_Request *mRequests;
        MPI_Status *status;
    };

  }
}

#endif // HEMELB_NET_H
