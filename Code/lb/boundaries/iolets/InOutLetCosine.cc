#include "lb/boundaries/iolets/InOutLetCosine.h"
#include "configuration/SimConfig.h"
#include "topology/NetworkTopology.h"

namespace hemelb
{
  namespace lb
  {
    namespace boundaries
    {
      namespace iolets
      {
        InOutLetCosine::InOutLetCosine() :
            InOutLet(), pressureMeanPhysical(0.0), pressureAmpPhysical(0.0), phase(0.0), period(1.0)
        {

        }

        void InOutLetCosine::DoIO(TiXmlElement *iParent, bool iIsLoading, configuration::SimConfig* iSimConfig)
        {
          iSimConfig->DoIOForCosineInOutlet(iParent, iIsLoading, this);
        }

        InOutLet* InOutLetCosine::Clone() const
        {
          InOutLetCosine* copy = new InOutLetCosine(*this);

          return copy;
        }

        InOutLetCosine::~InOutLetCosine()
        {

        }

        LatticeDensity InOutLetCosine::GetDensity(unsigned long time_step) const
        {
          double w = 2.0 * PI / period;
          return GetDensityMean() + GetDensityAmp() * cos(w * units->ConvertTimeStepToPhysicalUnits(time_step) + phase);
        }

        LatticeDensity InOutLetCosine::GetDensityMean() const
        {
          return units->ConvertPressureToLatticeUnits(pressureMeanPhysical) / Cs2;
        }

        LatticeDensity InOutLetCosine::GetDensityAmp() const
        {
          return units->ConvertPressureGradToLatticeUnits(pressureAmpPhysical) / Cs2;
        }

      }
    }
  }
}
