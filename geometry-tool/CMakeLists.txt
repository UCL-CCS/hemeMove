# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

cmake_minimum_required(VERSION 3.13)

project(HemeLbGmyTool)

# list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../CMake)
# include(Debug)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)

include(GNUInstallDirs)

find_package(ZLIB REQUIRED)
find_package(Boost REQUIRED COMPONENTS headers)
find_package(VTK 9 REQUIRED COMPONENTS
  CommonCore
  CommonDataModel
  CommonExecutionModel
  CommonMath
  FiltersCore
  FiltersGeneral
  FiltersSources
  IOXML
  WrappingPythonCore
)
find_package(CGAL REQUIRED)

if(SKBUILD)
  # Scikit-Build does not add your site-packages to the search path
  # automatically, so we need to add it _or_ the pybind11 specific
  # directory here.
  execute_process(
    COMMAND
      "${PYTHON_EXECUTABLE}" -c
      "import pybind11; print(pybind11.get_cmake_dir())"
    OUTPUT_VARIABLE _pybind11_cmake_dir
    OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND_ECHO STDOUT)
  list(APPEND CMAKE_PREFIX_PATH "${_pybind11_cmake_dir}")
endif()

find_package(pybind11 CONFIG REQUIRED)

set(HEMELB_CODE_DIR ${PROJECT_SOURCE_DIR}/../Code)
file(REAL_PATH "${HEMELB_CODE_DIR}" HEMELB_CODE_DIR)
set(model_dir HlbGmyTool/Model)
set(generation_dir HlbGmyTool/Model/Generation)

# First, our micro-HemeLB library (i.e. the bits we re-use).
set(hemelb_cpps
  util/Vector3D.cc
  geometry/SiteDataBare.cc
  io/writers/XdrFileWriter.cc
  io/writers/XdrMemWriter.cc
  io/writers/XdrWriter.cc
  io/writers/Writer.cc
  io/FILE.cc
)
list(TRANSFORM hemelb_cpps PREPEND "${HEMELB_CODE_DIR}/")
# While this is located in the extension source tree, in spirit it is
# part of hemelb::util::Vector3D's implementation.
list(APPEND hemelb_cpps "${generation_dir}/Index.cpp")

# This is a static archive as we want the objects to get pulled into
# our final single extension library to avoid the
# LD_LIBRARY_PATH/RPATH can of worms. But we require libuhlb to have
# the same POSITION_INDEPENDANT_CODE property as Generation, hence the
# function copy_picness below.
add_library(uhlb STATIC ${hemelb_cpps})
target_include_directories(uhlb PUBLIC ${HEMELB_CODE_DIR})

# Now the C++ implementation of the extension
set(gen_cpps
  Neighbours.cpp
  BuildCGALPolygon.cpp
  Block.cpp
  BlockWriter.cpp
  BufferPool.cpp
  GeometryGenerator.cpp
  GeometryWriter.cpp
  Domain.cpp
  Site.cpp
  InconsistentFluidnessError.cpp
  CylinderGenerator.cpp
  PolyDataGenerator.cpp
  SquareDuctGenerator.cpp
  Debug.cpp
  Wrap.cpp
)
list(TRANSFORM gen_cpps PREPEND "${generation_dir}/")

pybind11_add_module(Generation ${gen_cpps})

function(copy_picness src_tgt dest_tgt)
  set(pic POSITION_INDEPENDENT_CODE)
  get_target_property(picness ${src_tgt} ${pic})
  set_target_properties(${dest_tgt} PROPERTIES ${pic} ${picness})
endfunction()
copy_picness(Generation uhlb)

target_include_directories(Generation PRIVATE ${generation_dir})
target_link_libraries(Generation PRIVATE
  uhlb
  Boost::headers
  ZLIB::ZLIB
  Boost::headers
  CGAL::CGAL
  ${VTK_LIBRARIES}
)

install(TARGETS Generation DESTINATION ${model_dir})
