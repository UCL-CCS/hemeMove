# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

message(WARNING "The multiscale module is unmaintained - this may not work!")

find_hemelb_dependency(MPWide REQUIRED)

add_library(hemelb_multiscale OBJECT
  mpwide/MPWideIntercommunicator.cc
  )
target_link_libraries(hemelb_multiscale PUBLIC MPWide::MPWide)
