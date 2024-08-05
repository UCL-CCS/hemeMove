# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
add_test_lib(test-helpers
  FolderTestFixture.cc
  FourCubeBasedTestFixture.cc
  FourCubeLatticeData.cc
  HasCommsTestFixture.cc
  LabelledRequest.cc
  LatticeDataAccess.cc
  MockNetHelper.cc
  OneInOneOutSimConfig.cc
  RandomSourceTests.cc
  RecordingNet.cc
)
target_link_libraries(test-helpers PRIVATE tinyxml2::tinyxml2)
