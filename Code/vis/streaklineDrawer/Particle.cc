#include "Particle.h"

namespace hemelb
{
  namespace vis
  {
    namespace streaklinedrawer
    {
      Particle::Particle(float iX, float iY, float iZ, unsigned int iInletId) :
          x(iX), y(iY), z(iZ), vel(0), vx(0), vy(0), vz(0), inletID(iInletId)
      {
      }

      Particle::Particle(float iX, float iY, float iZ, float iVel, unsigned int iInletId) :
          x(iX), y(iY), z(iZ), vel(iVel), vx(0), vy(0), vz(0), inletID(iInletId)
      {
      }

    }
  }
}
