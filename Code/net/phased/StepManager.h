// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_NET_PHASED_STEPMANAGER_H
#define HEMELB_NET_PHASED_STEPMANAGER_H
#include <vector>
#include <set>
#include <map>
#include <string>
#include "reporting/Timers.h"
#include "net/IteratedAction.h"
#include "net/phased/Concern.h"
#include "net/phased/steps.h"

#include "log/Logger.h"
namespace hemelb
{
  namespace net
  {
    namespace phased
    {
      /***
       * The step manager provides ability to call, in sequence, methods on several
       * interested objects (Called "concerns") for each of several Steps,
       * in each of several Phases.
       * This provides a mechanism to allow concerns to schedule there activities
       * during and between each of the steps of an asynchronous communication pattern
       * There may be several sequences of asynchronous communication during a single
       * iteration of the SimulationController, by registering several phases, each with
       * it's own pre/post send-receive steps, etc.
       */
      class StepManager
      {
        public:

          typedef unsigned int Phase;
          typedef int MethodLabel;

          /***
           * Class abstracting methods which can be called on concerns
           * In principle, we would like to use function pointer classes to do this
           * But to keep things simpler, we are assuming indirection through a 'CallAction' method
           * on concerns, which handles dispatch to the appropriate method.
           */
          class Action
          {
            public:
              Concern * concern;
              MethodLabel method;
              std::string name;
              Action(Concern &concern, MethodLabel method) :
                  concern(&concern), method(method)
              {
              }
              Action(const Action & action) :
                  concern(action.concern), method(action.method)
              {
              }
              bool Call()
              {
                return concern->CallAction(method);
              }
          };

          typedef std::map<steps::Step, std::vector<Action> > Registry;

          /***
           * Construct a step manager
           * @param The number of phases, default 1.
           * @param timers, Record the times for the steps to this timers object, if given
           */
          StepManager(Phase phases = 1, reporting::Timers * timers = nullptr,
                      bool separate_concerns = false);

          /***
           * Register an action of a concern
           * @param phase Phase where an action should be called
           * @param step Step where an action should be called
           * @param concern Concern on which the method should be called
           * @param method label, indicating which method of the concern should be registered
           */
          void Register(Phase phase, steps::Step step, Concern & concern, MethodLabel method);

          /***
           * Register a concern for all of the steps typically used by an action,
           * syntactic sugar for the individual registrations
           * @param concern Concern which should be called, typically an IteratedActor
           * @param phase Phase for which this IteratedActor should be called.
           */
          void RegisterIteratedActorSteps(Concern &concern, Phase phase = 0);

          /***
           * Register a concern for all of the steps typically used to send/receive communications
           * syntactic sugar for the individual registrations (to Send, Receive, Wait)
           * @param concern Concern which should be called, typically a net::Net
           * @param phase Phase for which this concern should be used for comms
           */
          void RegisterCommsSteps(Concern &concern, Phase phase = 0);

          /***
           * Register a concern for all comms steps in all phases, typically the main net::Net
           * @param concern Concern to register, typically a net::Net
           */
          void RegisterCommsForAllPhases(Concern &concern);

          /***
           * Call the actions, concern by concern, for the given phase
           * @param phase
           */
          void CallActionsForPhaseSeparatedConcerns(Phase phase = 0);

          /***
           * Call the actions for the given phase
           * @param phase
           */
          void CallActionsForPhase(Phase phase = 0);

          /***
           * Call all registered actions for a given phase in a given step
           * @param step
           * @param phase (Default the first phase)
           */
          void CallActionsForStepForConcern(steps::Step step, Concern * concern, Phase phase = 0);

          /***
           * Call all registered actions for a given phase in a given step
           * @param step
           * @param phase (Default the first phase)
           */
          void CallActionsForStep(steps::Step step, Phase phase = 0);

          /***
           * Call all actions for all phases
           */
          void CallActions();

          /***
           * Call actions in a different order: do all steps for each concerns separately
           */
          void CallActionsSeparatedConcerns();

          /***
           * Call actions which do not belong to phases:
           * the special actions which occur before all phases and after all phases
           * @param step
           */
          void CallSpecialAction(steps::Step step);

          /***
           * Get the total number of registered concerns
           * @return total number of registered concerns
           */
          unsigned int ConcernCount() const;

          /***
           * Get the total number of registered actions
           * @return total number of registered actions
           */
          unsigned int ActionCount() const;

        private:
          std::vector<Registry> registry; // one registry for each phase
          std::vector<Concern*> concerns; // can't be a set as must be order-stable
          reporting::Timers *timers;
          void StartTimer(steps::Step step);
          void StopTimer(steps::Step step);
          const bool separate_concerns;

      };
    }
  }
}
#endif //ONCE
