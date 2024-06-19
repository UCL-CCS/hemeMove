# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

cmake_minimum_required(VERSION 3.13)

project(HemeLbGmyTool)

include(CTest)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../CMake)
include(Debug)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)

include(GNUInstallDirs)

find_package(Threads REQUIRED)
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
find_package(HDF5 REQUIRED COMPONENTS C)

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

set(CMAKE_CXX_VISIBILITY_PRESET default)
find_package(pybind11 CONFIG REQUIRED)

set(HEMELB_CODE_DIR ${PROJECT_SOURCE_DIR}/../Code)
file(REAL_PATH "${HEMELB_CODE_DIR}" HEMELB_CODE_DIR)
set(model_dir HlbGmyTool/Model)
set(gmy_dir ${model_dir}/GmyGeneration)
set(oct_dir ${model_dir}/OctGeneration)
set(common_dir ${model_dir}/CommonGeneration)

# First, our micro-HemeLB library (i.e. the bits we re-use).
set(hemelb_cpps
  util/Vector3D.cc
  geometry/SiteDataBare.cc
  io/formats/geometry.cc
  io/writers/XdrFileWriter.cc
  io/writers/XdrMemWriter.cc
  io/writers/XdrWriter.cc
  io/writers/Writer.cc
)
list(TRANSFORM hemelb_cpps PREPEND "${HEMELB_CODE_DIR}/")
# While this is located in the extension source tree, in spirit it is
# part of hemelb::util::Vector3D's implementation.
list(APPEND hemelb_cpps "${common_dir}/Index.cpp")

# This is a static archive as we want the objects to get pulled into
# our final single extension library to avoid the
# LD_LIBRARY_PATH/RPATH can of worms. But we require libuhlb to have
# the same POSITION_INDEPENDANT_CODE property as Generation, hence the
# function copy_picness below.
add_library(uhlb STATIC ${hemelb_cpps})
target_include_directories(uhlb PUBLIC ${HEMELB_CODE_DIR})
target_include_directories(uhlb PUBLIC ${common_dir})

# Now the C++ implementation of the extension
pybind11_add_module(CommonGeneration ${common_dir}/Wrap.cpp)
target_link_libraries(CommonGeneration PRIVATE uhlb)

# Older GMY file generator
set(gmy_gen_cpps
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
list(TRANSFORM gmy_gen_cpps PREPEND "${gmy_dir}/")

pybind11_add_module(GmyGeneration ${gmy_gen_cpps})

target_include_directories(GmyGeneration PRIVATE ${gmy_dir})
target_link_libraries(GmyGeneration PRIVATE
  uhlb
  Boost::headers
  ZLIB::ZLIB
  CGAL::CGAL
  vtksys
  vtkFiltersGeneral
  vtkCommonDataModel
  vtkIOXML
  vtkWrappingPythonCore
)

# Newer octree based one
set(oct_gen_cpps
  FloodFill.cpp
  H5.cpp
  MkCgalMesh.cpp
  PolyDataGenerator.cpp
  SectionTree.cpp
  SectionTreeBuilder.cpp
  SegmentFactory.cpp
  SurfaceVoxeliser.cpp
  TriTree.cpp
  TriangleSorter.cpp
)
list(TRANSFORM oct_gen_cpps PREPEND "${oct_dir}/")

add_library(octgmygen STATIC ${oct_gen_cpps})
target_include_directories(octgmygen PUBLIC ${oct_dir})
target_link_libraries(octgmygen
  PUBLIC
  uhlb
  PRIVATE
  vtksys
  vtkCommonDataModel
  Threads::Threads
  hdf5::hdf5
  Boost::headers
  ZLIB::ZLIB
  CGAL::CGAL
)

pybind11_add_module(OctGeneration ${oct_dir}/Wrap.cpp)

target_link_libraries(OctGeneration
  PRIVATE
  uhlb
  vtksys vtkCommonDataModel vtkWrappingPythonCore
  octgmygen
)

function(copy_picness src_tgt dest_tgt)
  set(pic POSITION_INDEPENDENT_CODE)
  get_target_property(picness ${src_tgt} ${pic})
  set_target_properties(${dest_tgt} PROPERTIES ${pic} ${picness})
endfunction()
copy_picness(GmyGeneration uhlb)
copy_picness(GmyGeneration octgmygen)

install(TARGETS CommonGeneration DESTINATION ${model_dir})
install(TARGETS GmyGeneration DESTINATION ${model_dir})
install(TARGETS OctGeneration DESTINATION ${model_dir})

if (BUILD_TESTING)
  find_package(Catch2 REQUIRED)
  include(Catch)

  set(test_cpps
    test.cpp
    TestResources/Meshes.cpp
    test_FloodFill.cpp
    test_H5.cpp
    test_MkCgalMesh.cpp
    test_Oct.cpp
    test_SectionTree.cpp
    test_SegmentFactory.cpp
    test_SurfaceVoxeliser.cpp
    test_TriangleSorter.cpp
  )
  list(TRANSFORM test_cpps PREPEND "${oct_dir}/")
  add_executable(hlbgmytool-octtest "${test_cpps}")
  target_link_libraries(hlbgmytool-octtest
    octgmygen
    vtksys
    vtkCommonDataModel
    vtkIOXML
    hdf5::hdf5
    Catch2::Catch2
  )

  catch_discover_tests(hlbgmytool-octtest
    PROPERTIES
    ENVIRONMENT "TESTRESOURCES=${PROJECT_SOURCE_DIR}/${oct_dir}/TestResources"
  )
endif()
