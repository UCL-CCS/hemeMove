# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
include_guard()

# This file allows you to set up a "pass" - a group of options that
# can be declared and then easily passed from one CMake project to
# another.

function(pass_group_declare GROUP)
  # Private: list of options
  SET(_${GROUP}_OPTIONS "" PARENT_SCOPE)
  # Private: list of the cached variables
  SET(_${GROUP}_CACHEVARS "" PARENT_SCOPE)
  # Private: list of variables to be forwarded from the super to the
  # code builds
  SET(_${GROUP}_FWDVARS "" PARENT_SCOPE)
  # Private: class with members for each option/variable
  SET(_${GROUP}_CLSDEFN "" PARENT_SCOPE)
endfunction()

#
# Declare an option (ON/OFF value)
#
function(pass_option GROUP NAME DESC DEFAULT)
  option(${NAME} ${DESC} ${DEFAULT})
  set(_${GROUP}_OPTIONS "${_${GROUP}_OPTIONS};${NAME}" PARENT_SCOPE)
  string(REPLACE "HEMELB_" "" BASENAME "${NAME}")
  set(_${GROUP}_CLSDEFN "${_${GROUP}_CLSDEFN}        static constexpr bool ${BASENAME} = \${${NAME}};
" PARENT_SCOPE)
endfunction()


# Types for the cachevar functions below are the allowed CMake types
# (BOOL, FILEPATH, PATH, STRING, INTERNAL) plus INT.
# These helpers convert from these to supported types in CMake/C++.

# Note that std::filesystem::path is not constexpr so use ct_string for those
function(_pass_type_to_cpp type out_type out_val)
  if (type STREQUAL "INT")
    set(${out_type} "int" PARENT_SCOPE)
    set(${out_val} "\${${NAME}}" PARENT_SCOPE)
  elseif (type STREQUAL "BOOL")
    set(${out_type} "bool" PARENT_SCOPE)
    set(${out_val} "\${${NAME}}" PARENT_SCOPE)
  else()
    set(${out_type} "ct_string" PARENT_SCOPE)
    set(${out_val} "\"\${${NAME}}\"" PARENT_SCOPE)
  endif()
endfunction()

# CMake doesn't do numbers, use string
function(_pass_type_to_cachevar type outvar)
  if (type STREQUAL "INT")
    set(${outvar} "STRING" PARENT_SCOPE)
  else()
    set(${outvar} ${type} PARENT_SCOPE)
  endif()
endfunction()

#
# Declare a cache variable
#
function(pass_cachevar GROUP NAME DEFAULT TYPE DESC)
  _pass_type_to_cachevar(${TYPE} cache_type)
  set(${NAME} ${DEFAULT} CACHE ${cache_type} ${DESC})
  set(_${GROUP}_CACHEVARS "${_${GROUP}_CACHEVARS};${NAME}" PARENT_SCOPE)
  string(REPLACE "HEMELB_" "" BASENAME "${NAME}")
  _pass_type_to_cpp(${TYPE} cpp_type cpp_val)
  set(_${GROUP}_CLSDEFN "${_${GROUP}_CLSDEFN}        static constexpr ${cpp_type} ${BASENAME} = ${cpp_val};
" PARENT_SCOPE)
endfunction()

function(pass_cachevar_choice GROUP NAME DEFAULT TYPE DESC)
  # ${ARGN} holds the choices
  _pass_type_to_cachevar(${TYPE} cache_type)
  set(${NAME} ${DEFAULT} CACHE ${cache_type} ${DESC})
  set(_${GROUP}_CACHEVARS "${_${GROUP}_CACHEVARS};${NAME}" PARENT_SCOPE)
  string(REPLACE "HEMELB_" "" BASENAME "${NAME}")
  _pass_type_to_cpp(${TYPE} cpp_type cpp_val)
  set(_${GROUP}_CLSDEFN "${_${GROUP}_CLSDEFN}        static constexpr ${cpp_type} ${BASENAME} = ${cpp_val};
" PARENT_SCOPE)
  set_property(CACHE ${NAME} PROPERTY STRINGS ${ARGN})
endfunction()


#
# Declare a variable to be forwarded from super build to code build
#
function(pass_var GROUP NAME)
  set(_${GROUP}_FWDVARS "${_${GROUP}_FWDVARS};${NAME}" PARENT_SCOPE)
endfunction()

#
# Function to produce the CMake defines to pass data from one build
# to an inner one. This will append to the output variable (so you
# may wish to set it to "" first).
#
function(pass_cmake_defines GROUP output)
  set(ans "${${output}}")
  foreach(OPT_NAME ${_${GROUP}_OPTIONS})
    list(APPEND ans "-D${OPT_NAME}=${${OPT_NAME}}")
  endforeach()

  foreach(V_NAME ${_${GROUP}_CACHEVARS})
    list(APPEND ans "-D${V_NAME}=${${V_NAME}}")
  endforeach()

  foreach(F_NAME ${_${GROUP}_FWDVARS})
    if (DEFINED ${F_NAME})
      list(APPEND ans "-D${F_NAME}=${${F_NAME}}")
    endif()
  endforeach()

  set(${output} ${ans} PARENT_SCOPE)
endfunction()

function(pass_get_defines GROUP output)
  string(CONFIGURE "${_${GROUP}_CLSDEFN}" tmp)
  SET(${output} "${tmp}" PARENT_SCOPE)
endfunction()
