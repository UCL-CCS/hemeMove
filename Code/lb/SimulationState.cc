// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "lb/SimulationState.h"

namespace hemelb::lb
{

    SimulationState::SimulationState(double timeStepLength, unsigned long endTS) :
        timeStepLength(timeStepLength), endTimeStep(endTS),
            isTerminating(false), stability(Stable)
    {
    }

    void SimulationState::Increment()
    {
      ++timeStep;
    }

    void SimulationState::Reset()
    {
      timeStep = startTimeStep;
    }

    void SimulationState::SetIsTerminating(bool value)
    {
      isTerminating = value;
    }

    void SimulationState::SetStability(Stability value)
    {
      stability = value;
    }

    bool SimulationState::IsTerminating() const
    {
      return isTerminating;
    }
    Stability SimulationState::GetStability() const
    {
      return stability;
    }

    void SimulationState::Report(reporting::Dict& dictionary)
    {
      dictionary.SetFormattedValue("TIME_STEP_LENGTH", "%lf", GetTimeStepLength());
      dictionary.SetIntValue("STEPS", GetTimeStepsElapsed());
      dictionary.SetIntValue("TOTAL_TIME_STEPS", endTimeStep - startTimeStep);

      switch (stability)
      {
        case lb::Unstable:
          dictionary.AddSectionDictionary("UNSTABLE");
          break;
        case lb::StableAndConverged:
          dictionary.AddSectionDictionary("SOLUTIONCONVERGED");
          break;
        default:
          break;
      }
    }
}
