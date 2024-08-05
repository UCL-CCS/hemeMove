// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "net/phased/StepManager.h"
#include <algorithm>

namespace hemelb
{
  namespace net
  {
    namespace phased
    {

      StepManager::StepManager(Phase phases, reporting::Timers *timers, bool separate_concerns) :
          registry(phases), concerns(), timers(timers), separate_concerns(separate_concerns)
      {
      }

      void StepManager::Register(Phase phase, steps::Step step, Concern & concern,
                                 MethodLabel method)
      {
        if (step == steps::BeginAll || step == steps::EndAll)
        {
          phase = 0; // special actions are always recorded in the phase zero registry
        }
        registry[phase][step].push_back(Action(concern, method));
        if (std::find(concerns.begin(), concerns.end(), &concern) == concerns.end())
        {
          concerns.push_back(&concern);
        }
      }

      void StepManager::RegisterIteratedActorSteps(Concern &concern, Phase phase)
      {
        for (int step = steps::BeginPhase; step <= steps::EndAll; step++)
        {
          if (step == steps::Receive || step == steps::Send || step == steps::Wait)
          {
            continue;
          }
          // C++ makes it completely annoying to iterate over values in an enum
          Register(phase, static_cast<steps::Step>(step), concern, step);
        }
      }

      void StepManager::RegisterCommsSteps(Concern &concern, Phase phase)
      {
        Register(phase, steps::Send, concern, steps::Send);
        Register(phase, steps::Receive, concern, steps::Receive);
        Register(phase, steps::Wait, concern, steps::Wait);
      }

      void StepManager::RegisterCommsForAllPhases(Concern &concern)
      {
        for (Phase phase = 0; phase < registry.size(); phase++)
        {
          RegisterCommsSteps(concern, phase);
        }
      }

      unsigned int StepManager::ConcernCount() const
      {
        return concerns.size();
      }

      unsigned int StepManager::ActionCount() const
      {
        unsigned int total = 0;
        for (std::vector<Registry>::const_iterator phase = registry.begin(); phase < registry.end();
            phase++)
        {
          for (Registry::const_iterator step = phase->begin(); step != phase->end(); step++)
          {
            for (std::vector<Action>::const_iterator action = step->second.begin();
                action != step->second.end(); action++)
            {
              total++;
            }
          }
        }
        return total;
      }

      void StepManager::CallActionsForPhase(Phase phase)
      {
        // It is assumed, that in the step enum, begin phase begins, and end phase ends, the steps which
        // must be called for a given phase.
        for (int step = steps::BeginPhase; step <= steps::EndPhase; step++)
        {
          CallActionsForStep(static_cast<steps::Step>(step), phase);
        }
      }

      void StepManager::CallActionsForPhaseSeparatedConcerns(Phase phase)
      {
        for (std::vector<Concern*>::iterator concern = concerns.begin(); concern != concerns.end();
            concern++)
        {
          for (int step = steps::BeginPhase; step <= steps::EndPhase; step++)
          {
            if (step == steps::Receive || step == steps::Send || step == steps::Wait)
            {
              // Call ALL comms actions for all concerns
              // Because, these concerns are net::Net objects, that do the actual send/receive/wait MPI calls
              /**
               * E.g:
               * if A is a concern, B is a concern, C is a concern
               * A is an it actor, B is an it actor,C is a net::Net
               * You want to go: A C A C A C B C B C B C
               */
              CallActionsForStep(static_cast<steps::Step>(step), phase);
            }
            else
            {
              // Call the actions only for THIS concern
              CallActionsForStepForConcern(static_cast<steps::Step>(step), *concern, phase);
            }
          }
        }
      }

      void StepManager::CallSpecialAction(steps::Step step)
      {
        // special actions are always recorded in the phase zero registry
        CallActionsForStep(static_cast<steps::Step>(step), 0);
      }

      void StepManager::CallActionsSeparatedConcerns()
      {
        CallSpecialAction(steps::BeginAll);
        for (Phase phase = 0; phase < registry.size(); phase++)
        {
          CallActionsForPhaseSeparatedConcerns(phase);
        }
        CallSpecialAction(steps::EndAll);
      }

      void StepManager::CallActions()
      {
        if (separate_concerns)
        {
          CallActionsSeparatedConcerns();
          return;
        }
        CallSpecialAction(steps::BeginAll);
        for (Phase phase = 0; phase < registry.size(); phase++)
        {
          CallActionsForPhase(phase);
        }
        CallSpecialAction(steps::EndAll);
      }

      void StepManager::CallActionsForStep(steps::Step step, Phase phase)
      {
        StartTimer(step);
        std::vector<Action> &actionsForStep = registry[phase][step];
        for (std::vector<Action>::iterator action = actionsForStep.begin();
            action != actionsForStep.end(); action++)
        {
          action->Call();
        }
        StopTimer(step);
      }

      void StepManager::CallActionsForStepForConcern(steps::Step step, Concern * concern,
                                                     Phase phase)
      {
        StartTimer(step);
        std::vector<Action> &actionsForStep = registry[phase][step];
        for (std::vector<Action>::iterator action = actionsForStep.begin();
            action != actionsForStep.end(); action++)
        {
          if (action->concern == concern)
          {
            action->Call();
          }
        }
        StopTimer(step);
      }

      void StepManager::StartTimer(steps::Step step)
      {
        if (!timers)
        {
          return;
        }
        if (step == steps::Wait)
        {
          timers->mpiWait().Start();
        }
        if (step == steps::Send)
        {
          timers->mpiSend().Start();
        }
      }

      void StepManager::StopTimer(steps::Step step)
      {
        if (!timers)
        {
          return;
        }
        if (step == steps::Wait)
        {
            timers->mpiWait().Stop();
        }
        if (step == steps::Send)
        {
            timers->mpiSend().Stop();
        }
      }
    }
  }
}
