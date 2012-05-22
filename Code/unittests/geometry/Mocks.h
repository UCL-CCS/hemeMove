#ifndef HEMELB_UNITTESTS_GEOMETRY_MOCKS_H
#define HEMELB_UNITTESTS_GEOMETRY_MOCKS_H
#include <vector>
#include <map>
#include <cstdlib>
#include <iostream>

#include "constants.h"
#include "mpiInclude.h"
#include "net/net.h"

namespace hemelb
{
  namespace unittests
  {
    namespace geometry
    {

      class NetMock
      {
        public:
          NetMock() :
              requiredReceipts(), requiredSends()
          {
          }
          void Receive()
          {
          }
          void Send()
          {
          }
          void Wait()
          {
          }
          void Dispatch()
          {
          }
          ~NetMock()
          {

          }

          class StoredCommunicationBase
          {
            public:
              virtual void RequireEnvelopeIdentity(const StoredCommunicationBase& other)=0;
              virtual void RequirePayloadIdentity(const StoredCommunicationBase& other)=0;
              virtual void Unpack(void *target)=0;
          };

          template<class T> class StoredCommunication : public StoredCommunicationBase
          {
            public:
              StoredCommunication(T* oPointer, unsigned int iCount, proc_t rank, const std::string &alabel = "") :
                  with(rank), payload(), count(iCount), label(alabel)
              {
                for (unsigned int i = 0; i < iCount; i++)
                {
                  payload.push_back(oPointer[i]); // by copy, slow, but ok cos this is just a mock
                }
              }
              void RequireEnvelopeIdentity(const StoredCommunicationBase& an_other)
              {
                const StoredCommunication<T> &other = dynamic_cast<const StoredCommunication<T> &>(an_other);
                std::stringstream message;
                message << "Destination rank unequal: " << label << " with message size " << count << std::endl;
                CPPUNIT_ASSERT_EQUAL_MESSAGE(message.str(), with, other.with);
                message.str("");
                message << "Message size unequal: " << label << " with core " << with << std::endl;
                CPPUNIT_ASSERT_EQUAL_MESSAGE(message.str(), count, other.count);
              }
              void RequirePayloadIdentity(const StoredCommunicationBase& an_other)
              {
                const StoredCommunication<T> &other = dynamic_cast<const StoredCommunication<T> &>(an_other);
                std::stringstream message;
                for (unsigned int i = 0; i < payload.size(); i++)
                {
                  message.str("");
                  message << "Element " << i << " unequal when comparing communicated data: " << label << std::flush;
                  CPPUNIT_ASSERT_EQUAL_MESSAGE(message.str(), payload[i], other.payload[i]);
                }
              }
              void Unpack(void *unknown_target)
              {
                T* target = static_cast<T*>(unknown_target);
                for (unsigned int i = 0; i < payload.size(); i++)
                {
                  target[i] = payload[i];
                }
              }
              proc_t with;
              std::vector<T> payload;
              unsigned int count;
              std::string label;
          };

          template<class T> void RequestSend(T* oPointer, int iCount, proc_t iToRank)
          {
            StoredCommunication<T> requirement(oPointer, iCount, iToRank);
            CPPUNIT_ASSERT(requiredSends.size() != 0);
            if (requiredSends.size() == 0)
            {
              return;
            }
            requiredSends.front()->RequireEnvelopeIdentity(requirement);
            requiredSends.front()->RequirePayloadIdentity(requirement);
            delete requiredSends.front();
            requiredSends.pop_front();
          }

          template<class T> void RequestReceive(T* oPointer, int iCount, proc_t iFromRank)
          {
            StoredCommunication<T> requirement(oPointer, iCount, iFromRank);
            CPPUNIT_ASSERT(requiredReceipts.size() != 0);
            if (requiredReceipts.size() == 0)
            {
              return;
            }
            requiredReceipts.front()->RequireEnvelopeIdentity(requirement);
            requiredReceipts.front()->Unpack(oPointer);
            delete requiredReceipts.front();
            requiredReceipts.pop_front();
          }

          /////////----------------------------------------------
          // Temporary copy in to make tests pass
          template<class T>
           void RequestSend(std::vector<T> &payload, proc_t toRank)
           {
             RequestSend(&payload[0], payload.size(), toRank);
           }

           template<class T>
           void RequestSend(T& value, proc_t toRank)
           {
             RequestSend(&value, 1, toRank);
           }

           template<class T>
           void RequestReceive(T& value, proc_t fromRank)
           {
             RequestReceive(&value, 1, fromRank);
           }

           template<class T>
           void RequestReceive(std::vector<T> &payload, proc_t toRank)
           {
             RequestReceive(&payload[0], payload.size(), toRank);
           }

           template<class T>
           void RequestGatherVReceive(std::vector<std::vector<T> > &buffer)
           {
             std::vector<int> displacements;
             std::vector<int> counts;
             for (typename std::vector<std::vector<T> >::iterator buffer_iterator = buffer.begin();
                 buffer_iterator != buffer.end(); buffer++)
             {
               displacements.push_back(&buffer_iterator->front() - &buffer.front());
               counts.push_back(buffer_iterator->size());
             }
             RequestGatherVReceive(&buffer.front(), &displacements.front(), &counts.front());
           }

           template<class T>
           void RequestGatherReceive(std::vector<T> &buffer)
           {

             RequestGatherReceive(&buffer.front());
           }

           template<class T>
           void RequestGatherSend(T& value, proc_t toRank)
           {
             RequestGatherSend(&value, toRank);
           }

           template<class T>
           void RequestGatherVSend(std::vector<T> &payload, proc_t toRank)
           {
             RequestGatherVSend(&payload.front(), payload.size(), toRank);
           }

           template<class T>
           void RequestGatherVSend(T* buffer, int count, proc_t toRank)
           {
           }

           template<class T>
           void RequestGatherReceive(T* buffer)
           {
           }

           template<class T>
           void RequestGatherSend(T* buffer, proc_t toRank)
           {
           }

           template<class T>
           void RequestGatherVReceive(T* buffer, int * displacements, int *counts)
           {
           }
           ///--------------------------------------------------

          template<class T> void RequireReceive(T* oPointer,
                                                unsigned int iCount,
                                                proc_t iFromRank,
                                                const std::string &label = "")
          {
            requiredReceipts.push_back(new StoredCommunication<T>(oPointer, iCount, iFromRank, label + " receive "));
          }
          template<class T> void RequireSend(T* oPointer,
                                             unsigned int iCount,
                                             proc_t iToRank,
                                             const std::string &label = "")
          {
            // Begin to fill up a buffer of Ts.
            // When RequestSends are given, they must match these, under operator ==, or an error will be raised.
            requiredSends.push_back(new StoredCommunication<T>(oPointer, iCount, iToRank, label + " send "));
          }
        private:
          std::deque<StoredCommunicationBase *> requiredReceipts;
          std::deque<StoredCommunicationBase *> requiredSends;
      };
    }
  }
}
#endif
