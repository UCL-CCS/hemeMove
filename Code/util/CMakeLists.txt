# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

add_library(hemelb_util OBJECT
  UnitConverter.cc clock.cc Vector3D.cc Vector3DHemeLb.cc Matrix3D.cc Bessel.cc)

if(LINUX_SCANDIR)
    target_compile_definitions(hemelb_util PRIVATE LINUX_SCANDIR)
endif()
