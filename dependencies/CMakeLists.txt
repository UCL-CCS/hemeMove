# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

cmake_minimum_required (VERSION 3.13)

project(HemeLBDependencies)

if(NOT DEFINED HEMELB_ROOT_DIR)
  set(HEMELB_ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/..)
endif()
if(NOT DEFINED HEMELB_DEPENDENCIES_PATH)
  set(HEMELB_DEPENDENCIES_PATH ${CMAKE_CURRENT_SOURCE_DIR})
endif()
if(NOT DEFINED HEMELB_DEPENDENCIES_INSTALL_PREFIX)
  set(HEMELB_DEPENDENCIES_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX} CACHE STRING "Path to install Dependencies")
endif()

list(INSERT CMAKE_MODULE_PATH 0 "${HEMELB_ROOT_DIR}/CMake")
list(INSERT CMAKE_MODULE_PATH 0 "${HEMELB_DEPENDENCIES_PATH}/Modules")
list(REMOVE_DUPLICATES CMAKE_MODULE_PATH)

list(INSERT CMAKE_PREFIX_PATH 0 "${HEMELB_DEPENDENCIES_INSTALL_PREFIX}")

include(ExternalProject)

include(GlobalOptions)
include(DependencyOptions)

find_package(MPI REQUIRED)

include(HemeLbDependency)

add_hemelb_dependency(Boost)
add_hemelb_dependency(tinyxml2)
add_hemelb_dependency(ParMETIS)
add_hemelb_dependency(CTemplate)
add_hemelb_dependency(ZLIB)

if (HEMELB_BUILD_TESTS)
  add_hemelb_dependency(Catch2)
endif()


if (HEMELB_BUILD_RBC)
  add_hemelb_dependency(HDF5)
  add_hemelb_dependency(VTK)
endif()

if (HEMELB_BUILD_MULTISCALE)
  add_hemelb_dependency(MPWide)
endif()

# Iff we are included from the superbuild, pass the list of active dependencies back.
get_directory_property(_have_parent PARENT_DIRECTORY)
if(_have_parent)
  set(HEMELB_PROJECT_DEPENDENCIES ${HEMELB_PROJECT_DEPENDENCIES} PARENT_SCOPE)
endif()
