// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#include "lb/iolets/InOutLetWomersleyVelocity.h"
#include "util/Bessel.h"

namespace hemelb::lb
{
      const InOutLetWomersleyVelocity::Complex InOutLetWomersleyVelocity::i = Complex(0, 1);
      const InOutLetWomersleyVelocity::Complex InOutLetWomersleyVelocity::iPowThreeHalves =
          pow(i, 1.5);

      InOutLet* InOutLetWomersleyVelocity::clone() const
      {
        InOutLet* copy = new InOutLetWomersleyVelocity(*this);
        return copy;
      }

      LatticeVelocity InOutLetWomersleyVelocity::GetVelocity(const LatticePosition& x,
                                                             const LatticeTimeStep t) const
      {
        LatticePosition displ = x - position;
        LatticeDistance z = Dot(displ, normal);
        Dimensionless r = sqrt(displ.GetMagnitudeSquared() - z * z);

        double omega = 2.0 * PI / period;
        LatticeDensity density = 1.0;

        Complex besselNumer = util::BesselJ0ComplexArgument(iPowThreeHalves * womersleyNumber * r
            / radius);
        Complex besselDenom = util::BesselJ0ComplexArgument(iPowThreeHalves * womersleyNumber);
        LatticeSpeed velocityMagnitude = std::real(pressureGradientAmplitude / (density * omega)
            * (1.0 - besselNumer / besselDenom) * exp(i * omega * double(t)));

        return normal * -velocityMagnitude;
      }

      const LatticePressureGradient& InOutLetWomersleyVelocity::GetPressureGradientAmplitude() const
      {
        return pressureGradientAmplitude;
      }

      void InOutLetWomersleyVelocity::SetPressureGradientAmplitude(
          const LatticePressureGradient& pressGradAmp)
      {
        pressureGradientAmplitude = pressGradAmp;
      }

      const LatticeTime& InOutLetWomersleyVelocity::GetPeriod() const
      {
        return period;
      }

      void InOutLetWomersleyVelocity::SetPeriod(const LatticeTime& per)
      {
        period = per;
      }

      const Dimensionless& InOutLetWomersleyVelocity::GetWomersleyNumber() const
      {
        return womersleyNumber;
      }

      void InOutLetWomersleyVelocity::SetWomersleyNumber(const Dimensionless& womNumber)
      {
        womersleyNumber = womNumber;
      }
    }
