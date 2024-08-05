// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_SIMULATIONSTATE_H
#define HEMELB_LB_SIMULATIONSTATE_H

#include "units.h"
#include "reporting/Reportable.h"
#include "units.h"

namespace hemelb::lb
{
    enum Stability
    {
      UndefinedStability = -1,
      Unstable = 0,
      Stable = 1,
      StableAndConverged = 2
    };

    class SimulationState : public reporting::Reportable
    {
    private:
        static constexpr LatticeTimeStep DEFAULT_INITIAL_TIME = 0;
        PhysicalTime timeStepLength;
        LatticeTimeStep timeStep = DEFAULT_INITIAL_TIME;
        LatticeTimeStep startTimeStep = DEFAULT_INITIAL_TIME;
        LatticeTimeStep endTimeStep;
        bool isTerminating;
        Stability stability;

        // To allow setting time stuff
        friend struct InitialConditionBase;

    public:
        SimulationState(double timeStepLength, unsigned long endTimeStep);
        ~SimulationState() noexcept override = default;

        void Increment();
        void Reset();
        void SetIsTerminating(bool value);
        void SetStability(Stability value);

        [[nodiscard]] inline LatticeTimeStep GetStartTimeStep() const {
            return startTimeStep;
        }
        [[nodiscard]] inline LatticeTimeStep GetEndTimeStep() const {
            return endTimeStep;
        }
        [[nodiscard]] inline LatticeTimeStep GetTimeStep() const {
            return timeStep;
        }
        [[nodiscard]] inline LatticeTimeStep GetTimeStepsElapsed() const {
            return timeStep - startTimeStep;
        }

        [[nodiscard]] bool IsTerminating() const;
        [[nodiscard]] Stability GetStability() const;

        [[nodiscard]] inline PhysicalTime GetTime() const
        {
          return GetTimeStepLength() * GetTimeStep();
        }
        [[nodiscard]] inline PhysicalTime GetTimeStepLength() const
        {
          return timeStepLength;
        }

        void Report(reporting::Dict& dictionary) override;
    };
}

#endif /* SIMULATIONSTATE_H_ */
