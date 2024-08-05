// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "net/IteratedAction.h"
#include "net/phased/StepManager.h"

namespace hemelb::net
{

    bool IteratedAction::CallAction(int action)
    {
      switch (static_cast<phased::steps::Step>(action))
      {
        case phased::steps::BeginPhase:
          RequestComms();
          return true;
        case phased::steps::PreSend:
          PreSend();
          return true;
        case phased::steps::PreWait:
          PreReceive();
          return true;
        case phased::steps::EndPhase:
          PostReceive();
          return true;
        case phased::steps::EndAll:
          EndIteration();
          return true;
        default:
          return false;
      }
    }

    void IteratedAction::RequestComms()
    {

    }

    void IteratedAction::PreSend()
    {

    }

    void IteratedAction::PreReceive()
    {

    }

    void IteratedAction::PostReceive()
    {

    }

    void IteratedAction::EndIteration()
    {

    }
}
