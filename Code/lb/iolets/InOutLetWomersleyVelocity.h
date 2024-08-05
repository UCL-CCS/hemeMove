// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_IOLETS_INOUTLETWOMERSLEYVELOCITY_H
#define HEMELB_LB_IOLETS_INOUTLETWOMERSLEYVELOCITY_H
#include "lb/iolets/InOutLetVelocity.h"
#include <complex>

namespace hemelb::lb
{
      /**
       * This class implements the InOutLetVelocity interface. It imposes a Womersley velocity profile
       * at a given inlet/outlet of a simulation of laminar flow in a cylinder of radius R and zero
       * average sinusoidal pressure gradient given by
       *
       *   pressureGradientAmplitude * sin(2*pi*t/period)
       *
       * where t is the current time. See physics validation paper for the analytical expression of
       * velocity implemented by GetVelocity(x, t) which is also a function of pressureGradientAmplitude,
       * radius, period, and womersleyNumber.
       *
       * If combined with a pressure iolet at the other end of the cylinder, it must be set to
       * zero pressure
       */
      class InOutLetWomersleyVelocity : public InOutLetVelocity
      {
        public:

          ~InOutLetWomersleyVelocity() override = default;

	  /**
           * Returns a copy of the current iolet. The caller is responsible for freeing that memory.
           *
           * @return copy of the current iolet
           */
          [[nodiscard]] InOutLet* clone() const override;

          /**
           * Get Womersley velocity for a given time and position.
           *
           * @param x lattice site position
           * @param t time
           * @return velocity
           */
          LatticeVelocity GetVelocity(const LatticePosition& x, const LatticeTimeStep t) const override;

          /**
           * Get the amplitude of the zero average pressure gradient sine wave imposed.
           *
           * @return pressure gradient amplitude
           */
          const LatticePressureGradient& GetPressureGradientAmplitude() const;

          /**
           * Set the amplitude of the zero average pressure gradient sine wave imposed.
           *
           * @param pressGradAmp pressure gradient amplitude
           */
          void SetPressureGradientAmplitude(const LatticePressureGradient& pressGradAmp);

          /**
           * Get the period of the zero average pressure gradient sine wave imposed
           *
           * @return sinusoidal pressure profile period
           */
          const LatticeTime& GetPeriod() const;

          /**
           * Set the period of the zero average pressure gradient sine wave imposed
           *
           * @param per sinusoidal pressure profile period
           */
          void SetPeriod(const LatticeTime& per);

          /**
           * Get the Womersley number characterising the pulsatile flow simulation
           *
           * @return Womersley number
           */
          const Dimensionless& GetWomersleyNumber() const;

          /**
           * Set the Womersley number characterising the pulsatile flow simulation
           *
           * @param womNumber Womersley number
           */
          void SetWomersleyNumber(const Dimensionless& womNumber);

        private:
          using Complex = std::complex<double>;
          static const Complex i;
          static const Complex iPowThreeHalves;
          LatticePressureGradient pressureGradientAmplitude; ///< See class documentation
          LatticeTime period; ///< See class documentation
          double womersleyNumber; ///< See class documentation
      };
}
#endif // HEMELB_LB_IOLETS_INOUTLETWOMERSLEYVELOCITY_H
