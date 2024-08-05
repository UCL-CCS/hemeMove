# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
include_guard()

include(PassOptions)
# This file defines all the CMake options and cache variables needed
# by HemeLB (main application). It can be included from either the the
# superbuild CMakeLists.txt or the app CMakeLists.txt

pass_group_declare(HEMELB)

#
# Specify the options
#
pass_option(HEMELB HEMELB_CALLGRIND "Compile unit tests that use callgrind" OFF)
pass_option(HEMELB HEMELB_USE_KRUEGER_ORDERING "Use Krueger's LB-IBM algorithm reodering" ON)
pass_option(HEMELB HEMELB_BUILD_DEBUGGER "Build the built in debugger" ON)
pass_option(HEMELB HEMELB_BUILD_COLLOIDS "Build the UNMAINTAINED colloids option" OFF)
# pass_option(HEMELB HEMELB_DEBUGGER_IMPLEMENTATION "Which implementation to use for the debugger" none)
# mark_as_advanced(HEMELB_DEBUGGER_IMPLEMENTATION)
pass_option(HEMELB HEMELB_VALIDATE_GEOMETRY "Validate geometry" OFF)
pass_option(HEMELB HEMELB_USE_ALL_WARNINGS_GNU "Show all compiler warnings on development builds (gnu-style-compilers)" ON)
pass_option(HEMELB HEMELB_DEPENDENCIES_SET_RPATH "Set runtime RPATH" ON)

if (CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")
  set(_default_sse3_flag ON)
else()
  set(_default_sse3_flag OFF)
endif()
pass_option(HEMELB HEMELB_USE_SSE3 "Use SSE3 intrinsics" ${_default_sse3_flag})
pass_option(HEMELB HEMELB_USE_VELOCITY_WEIGHTS_FILE "Use Velocity weights file" OFF)

pass_option(HEMELB HEMELB_SEPARATE_CONCERNS "Communicate for each concern separately" OFF)

if (HEMELB_BUILD_RBC)
  set(_default_kernel GuoForcingLBGK)
  set(_default_lattice D3Q19)
else()
  set(_default_kernel LBGK)
  set(_default_lattice D3Q15)
endif()

#
# Specify the variables
#
pass_cachevar(HEMELB HEMELB_EXECUTABLE "hemelb"
  STRING "File name of executable to produce")
pass_cachevar_choice(HEMELB HEMELB_LOG_LEVEL Info
  STRING "Log level"
  Critical Error Warning Info Debug Trace)

pass_cachevar_choice(HEMELB HEMELB_LATTICE ${_default_lattice}
  STRING "Select the lattice type to use"
  D3Q15 D3Q19 D3Q27 D3Q15i)
pass_cachevar_choice(HEMELB HEMELB_KERNEL ${_default_kernel}
  STRING "Select the kernel to use"
  LBGK EntropicAnsumali EntropicChik MRT TRT NNCY NNCYMOUSE NNC NNTPL GuoForcingLBGK)
pass_cachevar_choice(HEMELB HEMELB_WALL_BOUNDARY "SIMPLEBOUNCEBACK"
  STRING "Select the boundary conditions to be used at the walls"
  BFL GZS SIMPLEBOUNCEBACK JUNKYANG)
pass_cachevar_choice(HEMELB HEMELB_INLET_BOUNDARY "NASHZEROTHORDERPRESSUREIOLET"
  STRING "Select the boundary conditions to be used at the inlet"
  NASHZEROTHORDERPRESSUREIOLET LADDIOLET)
pass_cachevar_choice(HEMELB HEMELB_OUTLET_BOUNDARY "NASHZEROTHORDERPRESSUREIOLET"
  STRING "Select the boundary conditions to be used at the outlets"
  NASHZEROTHORDERPRESSUREIOLET LADDIOLET)
pass_cachevar_choice(HEMELB HEMELB_COMPUTE_ARCHITECTURE "AMDBULLDOZER"
  STRING "Select the architecture of the machine being used"
  INTELSANDYBRIDGE AMDBULLDOZER NEUTRAL ISBFILEVELOCITYINLET)
pass_cachevar(HEMELB HEMELB_POINTPOINT_IMPLEMENTATION Coalesce
  STRING "Point to point comms implementation, choose 'Coalesce', 'Separated', or 'Immediate'" )
pass_cachevar(HEMELB HEMELB_GATHERS_IMPLEMENTATION Separated
  STRING "Gather comms implementation, choose 'Separated', or 'ViaPointPoint'" )
pass_cachevar(HEMELB HEMELB_ALLTOALL_IMPLEMENTATION Separated
  STRING "Alltoall comms implementation, choose 'Separated', or 'ViaPointPoint'" )
pass_cachevar_choice(HEMELB HEMELB_STENCIL "FourPoint"
  STRING "HemeLB stencil type"
  TwoPoint ThreePoint FourPoint CosineApprox)

#
# Specify the variables requiring forwarding
#
pass_var(HEMELB CMAKE_INSTALL_PREFIX)
pass_var(HEMELB CMAKE_C_COMPILER)
pass_var(HEMELB CMAKE_CXX_COMPILER)

pass_var(HEMELB MPI_C_COMPILER)
pass_var(HEMELB MPI_CXX_COMPILER)
pass_var(HEMELB MPI_C_NO_INTERROGATE)
pass_var(HEMELB MPI_CXX_NO_INTERROGATE)
pass_var(HEMELB MPI_C_LIBRARIES)
pass_var(HEMELB MPI_C_INCLUDE_PATH)

pass_var(HEMELB BOOST_ROOT)
pass_var(HEMELB CTEMPLATE_USE_STATIC)
pass_var(HEMELB METIS_INCLUDE_DIR)
pass_var(HEMELB METIS_LIBRARY)
pass_var(HEMELB ParMETIS_INCLUDE_DIR)
pass_var(HEMELB ParMETIS_LIBRARY)
pass_var(HEMELB TINYXML2_INCLUDE_DIR)
pass_var(HEMELB TINYXML2_LIBRARIES)
