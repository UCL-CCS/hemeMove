// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UTIL_UNITCONVERTER_H
#define HEMELB_UTIL_UNITCONVERTER_H

#include "constants.h"
#include "units.h"
#include "util/Vector3D.h"
#include "util/Matrix3D.h"
#include "Exception.h"

namespace hemelb::util
{

    class UnitConverter
    {
    public:
        UnitConverter(PhysicalTime timeStep,
		      PhysicalDistance voxelSize, PhysicalPosition latticeOrigin,
		      PhysicalDensity fluidDensity, PhysicalPressure reference_pressure);

        LatticePressure ConvertPressureToLatticeUnits(PhysicalPressure pressure) const;
        LatticePressure ConvertPressureDifferenceToLatticeUnits(PhysicalStress pressure_diff) const;
        PhysicalPressure ConvertPressureToPhysicalUnits(LatticePressure pressure) const;

        LatticePressureGradient ConvertPressureGradientToLatticeUnits(PhysicalPressureGradient pg) const;
        PhysicalPressureGradient ConvertPressureGradientToPhysicalUnits(LatticePressureGradient pg) const;

        LatticeDistance ConvertDistanceToLatticeUnits(const PhysicalDistance& x) const;
        PhysicalDistance ConvertDistanceToPhysicalUnits(const LatticeDistance& x) const;

        LatticePosition ConvertPositionToLatticeUnits(const PhysicalPosition& x) const;
        PhysicalPosition ConvertPositionToPhysicalUnits(const LatticePosition& x) const;

        LatticeDisplacement ConvertDisplacementToLatticeUnits(const PhysicalDisplacement& x) const;
        PhysicalDisplacement ConvertDisplacementToPhysicalUnits(const LatticeDisplacement& x) const;

        LatticeSpeed ConvertSpeedToLatticeUnits(const PhysicalSpeed& v) const;
        PhysicalSpeed ConvertSpeedToPhysicalUnits(const LatticeSpeed& v) const;

        /**
         * Convert stress from physical to lattice units, using any rank of tensor
         */
        template<class InputType>
        InputType ConvertStressToLatticeUnits(InputType stress) const
        {
            using SCALAR = typename scalar_type<InputType>::type;
            return stress / SCALAR(latticePressure);
        }

        /**
         * Convert stress from lattice to physical units, using any rank of tensor
         */
        template<class InputType>
        InputType ConvertStressToPhysicalUnits(InputType stress) const
        {
            using SCALAR = typename scalar_type<InputType>::type;
            return stress * SCALAR(latticePressure);
        }

        /**
         * Convert a full stress tensor (including pressure and deviatoric components)
         * to physical units. Note how the diagonal is shifted by reference_pressure_Pa.
         *
         * @param stressTensor stress tensor in lattice units
         * @return stress tensor in physical units
         */
        Matrix3D ConvertFullStressTensorToPhysicalUnits(Matrix3D stressTensor) const
        {
          Matrix3D ret = stressTensor * latticePressure;
          ret.addDiagonal(reference_pressure_Pa);
          return ret;
        }

        /**
         * Convert a traction vector (force per unit area) to physical units. Note how a
         * reference_pressure_Pa*wallNormal component is added to account for the reference
         * pressure that was removed when converting the simulation input to lattice units.
         *
         * @param traction traction vector (computed the full stress tensor)
         * @param wallNormal wall normal at a given site
         * @return traction vector in physical units
         */
        template<class VectorType>
        Vector3D<VectorType> ConvertTractionToPhysicalUnits(
            Vector3D<VectorType> traction, const Vector3D<Dimensionless>& wallNormal) const
        {
            using SCALAR = typename scalar_type<VectorType>::type;
            Vector3D<VectorType> ret = traction * SCALAR(latticePressure);
            ret += wallNormal * reference_pressure_Pa;
            return ret;
        }

        /**
         * Templated to handle both absolute and directional velocity.
         * @param velocity
         * @return
         */
        template<class InputType>
        InputType ConvertVelocityToLatticeUnits(InputType velocity) const
        {
            using SCALAR = typename scalar_type<InputType>::type;
            return velocity / SCALAR(latticeSpeed);
        }

        /**
         * Templated to handle both absolute and directional velocity.
         * @param velocity
         * @return
         */
        template<class InputType>
        InputType ConvertVelocityToPhysicalUnits(InputType velocity) const
        {
            using SCALAR = typename scalar_type<InputType>::type;
            // convert velocity from lattice units to physical units (m/s)
            return velocity * SCALAR(latticeSpeed);
        }

        double ConvertPressureDifferenceToPhysicalUnits(distribn_t pressure_grad) const;

        LatticeTime ConvertTimeToLatticeUnits(const PhysicalTime& t) const;
        PhysicalTime ConvertTimeToPhysicalUnits(const LatticeTime& t) const;
        PhysicalTime ConvertTimeStepToPhysicalUnits(LatticeTimeStep time_step) const;

        /**
         * Converts a shear rate in lattice units into physical units
         * @param shearRate shear rate in lattice units (1/time_step_length)
         * @return shear rate in physical units (1/s)
         */
        PhysicalReciprocalTime ConvertShearRateToPhysicalUnits(
            LatticeReciprocalTime shearRate) const;

        inline const PhysicalDistance& GetVoxelSize() const
        {
          return latticeDistance;
        }

        inline const PhysicalTime& GetTimeStep() const {
            return latticeTime;
        }

        inline const PhysicalMass& GetMassScale() const {
            return latticeMass;
        }

        inline const PhysicalPressure& GetReferencePressure() const {
            return reference_pressure_Pa;
        }

        const PhysicalPosition& GetLatticeOrigin() const
        {
          return latticeOrigin;
        }
        LatticePosition GetPhysicalOrigin() const
        {
          return LatticePosition::Zero() - (latticeOrigin / latticeDistance);
        }

        template<typename T>
        T ConvertToLatticeUnits(std::string_view units, const T& value) const
        {
          double scale_factor;

          if (units == "m")
          {
            scale_factor = latticeDistance;
          }
          else if (units == "m/s/s")
          {
            scale_factor = latticeDistance / (latticeTime * latticeTime);
          }
          else if (units == "N")
          {
            // F = ma so Force = mass * length / time / time and Newton = Kg * metre / second / second
            scale_factor = latticeMass * latticeDistance / (latticeTime * latticeTime);
          }
          else if (units == "N/m")
          {
            scale_factor = latticeMass / (latticeTime * latticeTime);
          }
          else if (units == "Nm")
          {
            scale_factor = latticeMass * latticeDistance
              * latticeDistance / (latticeTime * latticeTime);
          }
          else if (units == "s")
          {
            scale_factor = latticeTime;
          }
          else if (units == "rad" || units == "dimensionless")
          {
            scale_factor = 1.;
          }
          else if (units == "m/s")
          {
            scale_factor = latticeDistance / latticeTime;
          }
          else if (units == "Pa")
          {
            scale_factor = latticeMass / (latticeDistance * latticeTime * latticeTime);
          }
          else if (units == "Pa/m")
          {
            scale_factor = latticeMass
                / (latticeDistance * latticeDistance * latticeTime * latticeTime);
          }
          else
          {
            throw Exception() << "Unknown units '" << units << "'";
          }
            using SCALAR = typename scalar_type<T>::type;
            return value / SCALAR(scale_factor);
        }

    private:
        PhysicalDistance latticeDistance; //!< Lattice displacement in physical units.
        PhysicalTime latticeTime;
        PhysicalMass latticeMass;
        PhysicalSpeed latticeSpeed; //!< Lattice displacement length divided by time step.
        PhysicalPosition latticeOrigin;
        PhysicalPressure latticePressure;
        PhysicalPressure reference_pressure_Pa;

        template <typename T>
        struct scalar_type {
            using type = T;
        };
        template <typename T>
        struct scalar_type<Vector3D<T>> {
            using type = T;
        };
    };

}

#endif /* HEMELB_UTIL_UNITCONVERTER_H */
