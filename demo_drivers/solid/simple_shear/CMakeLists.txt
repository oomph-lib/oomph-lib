# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered simple_shear subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(simple_shear C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

# Precompile the mesh as an OBJECT library
add_library(precompiled_mesh OBJECT precompiled_mesh.cc)
target_link_libraries(precompiled_mesh oomph::solid)
target_compile_features(precompiled_mesh PUBLIC cxx_std_17)

oomph_add_executable(
  NAME simple_shear
  SOURCES simple_shear.cc
  LIBRARIES oomph::solid oomph::constitutive oomph::meshes oomph::generic
            precompiled_mesh)

oomph_add_executable(
  NAME refineable_simple_shear
  SOURCES refineable_simple_shear.cc
  LIBRARIES oomph::solid oomph::constitutive oomph::meshes oomph::generic
            precompiled_mesh)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME solid.simple_shear
  TARGET_DEPENDENCIES simple_shear refineable_simple_shear
  EXTRA_REQUIRES
  LABELS solid simple_shear)

# This test needs more than the default allotted time to complete
set_tests_properties(solid.simple_shear PROPERTIES TIMEOUT 10000)

# TODO: Needs to be patched
set_tests_properties(solid.simple_shear PROPERTIES DISABLED YES)

message(VERBOSE "Leaving simple_shear subdirectory")
# ------------------------------------------------------------------------------
