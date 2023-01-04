# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered octree_test subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(octree_test C CXX Fortran)
if(NOT OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME octree_test
  SOURCES octree_test.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME poisson.octree_test
  TARGET_DEPENDENCIES octree_test
  EXTRA_REQUIRES
  LABELS octree poisson)

message(VERBOSE "Leaving octree_test subdirectory")
# ------------------------------------------------------------------------------
