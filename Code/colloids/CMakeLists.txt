# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
message(FATAL_ERROR "The colloids module is unmaintained - see README")
add_library(hemelb_colloids OBJECT
  BoundaryConditions.cc
  BodyForces.cc
  ParticleMpiDatatypes.cc
  ParticleSet.cc
  Particle.cc
  PersistedParticle.cc
  ColloidController.cc
  )
target_link_libraries(hemelb_colloids
  hemelb_io
  hemelb_configuration
  hemelb_lb
  )
target_compile_definitions(hemelb_colloids PUBLIC HEMELB_BUILD_COLLOIDS)
