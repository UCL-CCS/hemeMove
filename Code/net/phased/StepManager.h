#ifndef HEMELB_NET_PHASED_STEPMANAGER_H
#define HEMELB_NET_PHASED_STEPMANAGER_H
#include <vector>
#include <map>
#include <set>
#include "net/IteratedAction.h"
#include "net/phased/Concern.h"
#include "net/phased/steps.h"
namespace hemelb
{
  namespace net
  {
    namespace phased
    {


      class StepManager
      {
        public:

          typedef unsigned int Phase;
          typedef int MethodLabel;

          class Action
          {
            public:
              Concern & concern;
              MethodLabel method;
              Action(Concern &concern, MethodLabel method) :
                  concern(concern), method(method)
              {
              }
          };

          typedef std::pair<steps::Step, Action> BoundAction;
          typedef std::multimap<steps::Step, Action> Registry;

          StepManager(Phase phases = 1);

          void Register(Phase phase, steps::Step step, Concern & concern, MethodLabel method);

          void RegisterIteratedActorSteps(Concern &concern, Phase phase = 0);

          void RegisterCommsSteps(Concern &concern, Phase phase = 0);

          unsigned int ConcernCount() const;

          unsigned int ActionCount() const;

        private:
          std::vector<Registry> registry;

      };
    }
  }
}
#endif //ONCE
