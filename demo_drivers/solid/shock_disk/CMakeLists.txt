# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered shock_disk subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(shock_disk C CXX Fortran)
if(NOT OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME shock_disk
  SOURCES shock_disk.cc
  LIBRARIES oomph::solid oomph::constitutive oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME solid.shock_disk
  TARGET_DEPENDENCIES shock_disk
  EXTRA_REQUIRES
  LABELS solid shock_disk)

message(VERBOSE "Leaving shock_disk subdirectory")
# ------------------------------------------------------------------------------
