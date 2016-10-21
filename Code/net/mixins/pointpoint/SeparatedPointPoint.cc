
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "net/mixins/pointpoint/SeparatedPointPoint.h"
#include "log/Logger.h"

namespace hemelb
{
  namespace net
  {

    void SeparatedPointPoint::EnsureEnoughRequests(size_t count)
    {
      if (requests->size() < count)
      {
        requests->resize(count);
      }
    }

    void SeparatedPointPoint::ReceivePointToPoint()
    {
      EnsurePreparedToSendReceive();
      proc_t m = 0;

      for (std::map<proc_t, ProcComms>::iterator it = receiveProcessorComms.begin(); it != receiveProcessorComms.end();
          ++it)
      {
        for (ProcComms::iterator request = it->second.begin(); request != it->second.end(); request++)
        {
	  requests->set(m,
			communicator->IrecvImpl(request->Pointer,
						request->Count,
						request->Type,
						it->first,
						10));
          ++m;
        }
      }

    }

    // Makes sure the MPI_Datatypes for sending and receiving have been created for every neighbour.
    void SeparatedPointPoint::EnsurePreparedToSendReceive()
    {
      if (sendReceivePrepped)
      {
        return;
      }

      count_sends=0;
      count_receives=0;
      for (std::map<proc_t, ProcComms>::iterator it = sendProcessorComms.begin(); it != sendProcessorComms.end(); ++it)
      {
        count_sends += it->second.size();
      }

      for (std::map<proc_t, ProcComms>::iterator it = receiveProcessorComms.begin(); it != receiveProcessorComms.end();
          ++it)
      {
        count_receives += it->second.size();
      }

      EnsureEnoughRequests(count_sends+count_receives);

      sendReceivePrepped = true;
    }

    void SeparatedPointPoint::SendPointToPoint()
    {

      // Make sure the MPI datatypes have been created.
      EnsurePreparedToSendReceive();
      proc_t m = 0;

      for (std::map<proc_t, ProcComms>::iterator it = sendProcessorComms.begin(); it != sendProcessorComms.end(); ++it)
      {
        for (ProcComms::iterator request = it->second.begin(); request != it->second.end(); request++)
        {
	  requests->set(count_receives+m,
			communicator->IsendImpl(request->Pointer,
						request->Count,
						request->Type,
						it->first,
						10));
          ++m;
        }
      }
    }

    /*!
     Free the allocated data.
     */
    SeparatedPointPoint::~SeparatedPointPoint()
    {
    }

    void SeparatedPointPoint::WaitPointToPoint()
    {
      requests->WaitAll();

      receiveProcessorComms.clear();
      sendProcessorComms.clear();
      sendReceivePrepped = false;

    }
  }
}
