# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

add_library(hemelb_log OBJECT Logger.cc)

if(HAVE_RUSAGE)
    target_compile_definitions(hemelb_log PRIVATE HAVE_RUSAGE)
endif()
