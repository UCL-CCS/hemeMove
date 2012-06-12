#ifndef HEMELB_NET_PHASED_STEPS_H
#define HEMELB_NET_PHASED_STEPS_H

namespace hemelb
{
  namespace net
  {
    namespace phased
    {
      namespace steps
      {
        enum Step
        {
          BeginAll = -1, // Called only before first phase
          BeginPhase = 0,
          Receive = 1,
          PreSend = 2,
          Send = 3,
          PreWait = 4,
          Wait = 5,
          EndPhase = 6,
          EndAll = 7 // Called only after final phase
        };
      }
    }
  }
}
#endif
