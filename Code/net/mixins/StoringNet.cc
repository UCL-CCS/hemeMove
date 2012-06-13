#include "net/mixins/StoringNet.h"
namespace hemelb
{
  namespace net
  {

    void StoringNet::RequestSendImpl(void* pointer, int count, proc_t rank, MPI_Datatype type)
    {
      if (count > 0)
      {
        sendProcessorComms[rank].push_back(SimpleRequest(pointer, count, type, rank));
      }
    }

    void StoringNet::RequestReceiveImpl(void* pointer, int count, proc_t rank, MPI_Datatype type)
    {
      if (count > 0)
      {
        receiveProcessorComms[rank].push_back(SimpleRequest(pointer, count, type, rank));
      }
    }

    /*
     * Blocking gathers are implemented in MPI as a single call for both send/receive
     * But, here we separate send and receive parts, since this interface may one day be used for
     * nonblocking collectives.
     */

    void StoringNet::RequestGatherVSendImpl(void* buffer, int count, proc_t toRank, MPI_Datatype type)
    {
      gatherVSendProcessorComms[toRank].push_back(SimpleRequest(buffer, count, type, toRank));
    }

    void StoringNet::RequestGatherReceiveImpl(void* buffer, MPI_Datatype type)
    {
      /*
       * Dummy rank to ScalarRequest of zero -- gathers always receive to the core where the receive request is made.
       * Avoids defining another type of request.
       */
      gatherReceiveProcessorComms.push_back(ScalarRequest(buffer, type, 0));
    }

    void StoringNet::RequestGatherSendImpl(void* buffer, proc_t toRank, MPI_Datatype type)
    {
      gatherSendProcessorComms[toRank].push_back(ScalarRequest(buffer, type, toRank));
    }

    void StoringNet::RequestGatherVReceiveImpl(void* buffer, int * displacements, int *counts, MPI_Datatype type)
    {
      gatherVReceiveProcessorComms.push_back(GatherVReceiveRequest(buffer, displacements, counts, type));
    }
  }
}
