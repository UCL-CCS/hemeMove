// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "util/UnitConverter.h"
#include "constants.h"

namespace hemelb::util
{

    UnitConverter::UnitConverter(PhysicalTime timeStep,
				 PhysicalDistance voxelSize, PhysicalPosition latticeOrigin,
				 PhysicalDensity fluidDensity, PhysicalPressure reference_pressure) :
        latticeDistance(voxelSize), latticeTime(timeStep),
            latticeMass(fluidDensity * voxelSize * voxelSize * voxelSize),
            latticeSpeed(voxelSize / latticeTime), latticeOrigin(latticeOrigin),
            latticePressure(latticeMass / (latticeDistance * latticeTime * latticeTime)),
            reference_pressure_Pa(reference_pressure)
    {

    }

    LatticePressure UnitConverter::ConvertPressureToLatticeUnits(PhysicalPressure pressure) const
    {
      return Cs2 + ConvertPressureDifferenceToLatticeUnits(pressure - reference_pressure_Pa);
    }

    PhysicalPressure UnitConverter::ConvertPressureToPhysicalUnits(LatticePressure pressure) const
    {
      return reference_pressure_Pa + ConvertPressureDifferenceToPhysicalUnits(pressure - Cs2);
    }

    LatticePressure UnitConverter::ConvertPressureDifferenceToLatticeUnits(
        PhysicalPressure pressure_diff) const
    {
      return pressure_diff / latticePressure;
    }

    PhysicalPressure UnitConverter::ConvertPressureDifferenceToPhysicalUnits(
        LatticePressure pressure_diff) const
    {
      return pressure_diff * latticePressure;
    }

    LatticePressureGradient UnitConverter::ConvertPressureGradientToLatticeUnits(PhysicalPressureGradient pg) const {
        return pg * latticeDistance / latticePressure;
    }

    PhysicalPressureGradient UnitConverter::ConvertPressureGradientToPhysicalUnits(LatticePressureGradient pg) const {
        return pg * latticePressure / latticeDistance;
    }

    PhysicalReciprocalTime UnitConverter::ConvertShearRateToPhysicalUnits(
        LatticeReciprocalTime shearRate) const
    {
      return shearRate / latticeTime;
    }
    LatticeDistance UnitConverter::ConvertDistanceToLatticeUnits(const PhysicalDistance& x) const
    {
      return x / latticeDistance;
    }
    PhysicalDistance UnitConverter::ConvertDistanceToPhysicalUnits(const LatticeDistance& x) const
    {
      return x * latticeDistance;
    }
    LatticePosition UnitConverter::ConvertPositionToLatticeUnits(const PhysicalPosition& x) const
    {
      return (x - latticeOrigin) / latticeDistance;
    }
    PhysicalPosition UnitConverter::ConvertPositionToPhysicalUnits(const LatticePosition& x) const
    {
      return x * latticeDistance + latticeOrigin;
    }
    LatticeDisplacement UnitConverter::ConvertDisplacementToLatticeUnits(const PhysicalDisplacement& x) const
    {
        return x / latticeDistance;
    }
    PhysicalDisplacement UnitConverter::ConvertDisplacementToPhysicalUnits(const LatticeDisplacement& x) const
    {
        return x * latticeDistance;
    }
    LatticeSpeed UnitConverter::ConvertSpeedToLatticeUnits(const PhysicalSpeed& v) const
    {
      return v / latticeSpeed;
    }
    PhysicalSpeed UnitConverter::ConvertSpeedToPhysicalUnits(const LatticeSpeed& v) const
    {
      return v * latticeSpeed;
    }
    LatticeTime UnitConverter::ConvertTimeToLatticeUnits(const PhysicalTime& t) const
    {
      return t / latticeTime;
    }
    PhysicalTime UnitConverter::ConvertTimeToPhysicalUnits(const LatticeTime& t) const
    {
      return t * latticeTime;
    }
    PhysicalTime UnitConverter::ConvertTimeStepToPhysicalUnits(LatticeTimeStep time_step) const
    {
      return LatticeTime(time_step) * latticeTime;
    }

  }
