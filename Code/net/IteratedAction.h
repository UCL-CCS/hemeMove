// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_NET_ITERATEDACTION_H
#define HEMELB_NET_ITERATEDACTION_H
#include "net/phased/Concern.h"
namespace hemelb::net
  {
    class IteratedAction : public phased::Concern
    {
      public:
        /***
         * When an iterated actor is called through the phased::StepManager mechanism,
         * this method dispatches to the individual step methods
         * @param action Enumeration indicating the step
         * @return True if an action was successfully called for the step
         */
        bool CallAction(int action) final;
        //~IteratedAction() override = default;
        virtual void RequestComms();
        virtual void PreSend();
        virtual void PreReceive();
        virtual void PostReceive();
        virtual void EndIteration();
    };
  }

#endif /* HEMELB_NET_ITERATEDACTION_H */
