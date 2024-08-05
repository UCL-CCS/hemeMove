# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
include_guard()

find_file(TINYXML2_TARBALL tinyxml2-9.0.0.tar.gz 
  DOC "Path to download tinyxml2 (can be url http://)"
  PATHS ${HEMELB_DEPENDENCIES_PATH}/distributions
  )
if(NOT TINYXML_TARBALL)
  message("No tinyxml2 source found, will download.")
  set(TINYXML2_TARBALL https://github.com/leethomason/tinyxml2/archive/refs/tags/9.0.0.tar.gz
    CACHE STRING "Path to download TinyXML2 (can be local file://)" FORCE)
endif()
ExternalProject_Add(
  dep_tinyxml2
  INSTALL_DIR ${HEMELB_DEPENDENCIES_INSTALL_PREFIX}
  URL ${TINYXML2_TARBALL}
  CONFIGURE_COMMAND cmake <SOURCE_DIR>
  -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
  -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
  -Dtinyxml2_BUILD_TESTING=OFF
  -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
  -DCMAKE_CXX_FLAGS_RELEASE=${CMAKE_CXX_FLAGS_RELEASE}
  BUILD_COMMAND make
)
