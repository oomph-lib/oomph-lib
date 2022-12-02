# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered circular_disk subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(circular_disk C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME circular_disk
  SOURCES circular_disk.cc
  LIBRARIES oomph::foeppl_von_karman oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME foeppl_von_karman.circular_disk
  TARGET_DEPENDENCIES circular_disk
  EXTRA_REQUIRES
  LABELS foeppl_von_karman circular_disk)

message(VERBOSE "Leaving circular_disk subdirectory")
# ------------------------------------------------------------------------------
