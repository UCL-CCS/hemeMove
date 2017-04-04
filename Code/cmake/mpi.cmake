
# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
# ------MPI------------------
# Require MPI for this project:
find_package(MPI REQUIRED)
set(CMAKE_CXX_COMPILE_FLAGS "${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS}")
set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} ${CMAKE_CXX_LINK_FLAGS}")
include_directories(${MPI_INCLUDE_PATH})

function(TEST_MPI_VERSION_EQUAL ver output_var)
  set(CMAKE_REQUIRED_FLAGS "${CXX11_FLAGS} ${CMAKE_CXX_COMPILE_FLAGS}")
  set(CMAKE_REQUIRED_INCLUDES "${MPI_INCLUDE_PATH}")
  set(CMAKE_REQUIRED_LIBRARIES "${MPI_CXX_LIBRARIES}")
  CHECK_CXX_SOURCE_COMPILES("#include <mpi.h>
int main(int argc, char* argv[]) {
  static_assert(MPI_VERSION == ${ver}, \"\");
}" HAVE_MPI_STANDARD_VERSION_${ver})
  SET(${output_var} ${HAVE_MPI_STANDARD_VERSION_${ver}} PARENT_SCOPE)
endfunction()

function(TEST_MPI_SUBVERSION_EQUAL ver output_var)
  set(CMAKE_REQUIRED_FLAGS "${CXX11_FLAGS} ${CMAKE_CXX_COMPILE_FLAGS}")
  set(CMAKE_REQUIRED_INCLUDES "${MPI_INCLUDE_PATH}")
  set(CMAKE_REQUIRED_LIBRARIES "${MPI_CXX_LIBRARIES}")
  CHECK_CXX_SOURCE_COMPILES("#include <mpi.h>
int main(int argc, char* argv[]) {
  static_assert(MPI_SUBVERSION == ${ver}, \"\");
}" HAVE_MPI_STANDARD_SUBVERSION_${ver})
  SET(${output_var} ${HAVE_MPI_STANDARD_SUBVERSION_${ver}} PARENT_SCOPE)
endfunction()

function(GET_MPI_VERSION output_var)
  # Future proof as MPI 4.0 will appear eventually
  foreach(version RANGE 1 4)
    TEST_MPI_VERSION_EQUAL(${version} tmp)
    if(${tmp})
      SET(${output_var} ${version} PARENT_SCOPE)
      return()
    endif()
  endforeach()
  message(FATAL_ERROR "Could not determine MPI_VERSION")
endfunction()

function(GET_MPI_SUBVERSION output_var)
  foreach(version RANGE 9)
    TEST_MPI_SUBVERSION_EQUAL(${version} tmp)
    if(${tmp})
      SET(${output_var} ${version} PARENT_SCOPE)
      return()
    endif()
  endforeach()
  message(FATAL_ERROR "Could not determine MPI_SUBVERSION")
endfunction()

# CMake doens't seem to properly pass through the C++11 flag
# This is probably wrong
set(CMAKE_REQUIRED_FLAGS ${CXX11_FLAGS})
set(CMAKE_REQUIRED_QUIET 1)

# Put this in here as this is the key feature needed by the test programs
CHECK_CXX_SOURCE_COMPILES("int main(int, char**) {
static_assert(true, \"must not fail\");
}"
  HAVE_STATIC_ASSERT)

if(NOT HAVE_STATIC_ASSERT)
  message(FATAL_ERROR "No static_assert!")
endif()

GET_MPI_VERSION(MPI_STANDARD_VERSION_MAJOR)
GET_MPI_SUBVERSION(MPI_STANDARD_VERSION_MINOR)
SET(MPI_STANDARD_VERSION "${MPI_STANDARD_VERSION_MAJOR}.${MPI_STANDARD_VERSION_MINOR}")
message(STATUS "MPI library claims to implement standard version ${MPI_STANDARD_VERSION}")

unset(CMAKE_REQUIRED_FLAGS)
unset(CMAKE_REQUIRED_QUIET)
