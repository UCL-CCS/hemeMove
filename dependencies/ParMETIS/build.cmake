# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
include_guard()

find_file(ParMETIS_TARBALL parmetis-4.0.2.tar.gz 
  DOC "Path to download ParMETIS (can be url http://)"
  PATHS ${HEMELB_DEPENDENCIES_PATH}/distributions
  )
if(NOT ParMETIS_TARBALL)
  message("No ParMETIS source found, will download.")
  set(ParMETIS_TARBALL https://www.cp2k.org/static/downloads/parmetis-4.0.2.tar.gz 
    CACHE STRING "Path to download ParMETIS (can be local file://)" FORCE)
endif()
  
set(ParMETIS_CC ${MPI_C_COMPILER} CACHE STRING "MPI Compiler to use for ParMETIS, leave blank to let ParMETIS guess")
set(ParMETIS_CXX ${MPI_CXX_COMPILER} CACHE STRING "MPI Compiler to use for ParMETIS, leave blank to let ParMETIS guess")
if(ParMETIS_CC)
  set(ParMETIS_CC_OPTION cc=${ParMETIS_CC})
endif()
if(ParMETIS_CXX)
  set(ParMETIS_CXX_OPTION cxx=${ParMETIS_CXX})
endif()
ExternalProject_Add(
  dep_ParMETIS
  INSTALL_DIR ${HEMELB_DEPENDENCIES_INSTALL_PREFIX}
  URL ${ParMETIS_TARBALL}
  CONFIGURE_COMMAND make config prefix=<INSTALL_DIR> ${ParMETIS_CC_OPTION} ${ParMETIS_CXX_OPTION} && 
  cd metis && 
  make config prefix=<INSTALL_DIR> ${ParMETIS_CC_OPTION} ${ParMETIS_CXX_OPTION}
  BUILD_IN_SOURCE 1
  BUILD_COMMAND make -j${HEMELB_SUBPROJECT_MAKE_JOBS} MAKEFLAGS= #Workaround for bug in ParMETIS makefile treating switches as targets
  INSTALL_COMMAND make install MAKEFLAGS= && cd metis && make install MAKEFLAGS=
)
