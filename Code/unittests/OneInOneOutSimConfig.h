#ifndef HEMELB_UNITTESTS_ONEINONEOUTSIMCONFIG_H
#define HEMELB_UNITTESTS_ONEINONEOUTSIMCONFIG_H

#include "configuration/SimConfig.h"

namespace hemelb
{
  namespace unittests
  {
    class OneInOneOutSimConfig : public configuration::SimConfig
    {
      public:
        OneInOneOutSimConfig() :
          configuration::SimConfig()
        {
          lb::boundaries::iolets::InOutLetCosine* inlet = new lb::boundaries::iolets::InOutLetCosine();
          inlet->PressureAmpPhysical = 1.0;
          inlet->PressureMeanPhysical = 80.0;
          inlet->Phase = PI;

          Inlets.push_back(inlet);

          lb::boundaries::iolets::InOutLetCosine* outlet = new lb::boundaries::iolets::InOutLetCosine();
          outlet->PressureAmpPhysical = 0.0;
          outlet->PressureMeanPhysical = 80.0;
          outlet->Phase = 0.0;

          Outlets.push_back(outlet);

          TotalTimeSteps = 10000;
          TimeStepLength = 60.0/(70.0*1000.0);
        }
    };
  }
}

#endif /* HEMELB_UNITTESTS_ONEINONEOUTSIMCONFIG_H */
