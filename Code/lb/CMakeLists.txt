# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

add_library(hemelb_lb OBJECT
  iolets/BoundaryCommunicator.cc iolets/BoundaryComms.cc iolets/BoundaryValues.cc
  iolets/InOutLet.cc
  iolets/InOutLetCosine.cc iolets/InOutLetFile.cc
  iolets/InOutLetMultiscale.cc
  iolets/InOutLetVelocity.cc
  iolets/InOutLetParabolicVelocity.cc iolets/InOutLetWomersleyVelocity.cc iolets/InOutLetFileVelocity.cc
  IncompressibilityChecker.cc
        kernels/DHumieresD3Q15MRTBasis.cc kernels/DHumieresD3Q19MRTBasis.cc
  kernels/AbstractRheologyModel.cc kernels/CarreauYasudaRheologyModel.cc
  kernels/CassonRheologyModel.cc kernels/TruncatedPowerLawRheologyModel.cc
  MacroscopicPropertyCache.cc SimulationState.cc StabilityTester.cc
  InitialCondition.cc
        lb.cc
  )
