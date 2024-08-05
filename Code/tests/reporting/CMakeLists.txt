# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
add_test_lib(test_reporting
  ReporterTests.cc
  TimerTests.cc
)
target_link_libraries(test_reporting PRIVATE CTemplate::CTemplate)
