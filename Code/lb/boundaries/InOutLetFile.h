#ifndef HEMELB_LB_BOUNDARIES_INOUTLETFILE_H
#define HEMELB_LB_BOUNDARIES_INOUTLETFILE_H

#include "lb/boundaries/InOutLetCycle.h"

namespace hemelb
{
  namespace lb
  {
    namespace boundaries
    {

      class InOutLetFile : public InOutLetCycle<0, false>
      {
        public:
          InOutLetFile();
          virtual ~InOutLetFile();

          virtual void DoIO(TiXmlElement *iParent, bool iIsLoading, SimConfig* iSimConfig);
          virtual InOutLet* Clone();

          virtual void CalculateCycle(std::vector<distribn_t> &density_cycle,
                                      const SimulationState *iState);

          std::string PressureFilePath;
      };

    }
  }
}

#endif /* HEMELB_LB_BOUNDARIES_INOUTLETFILE_H */
