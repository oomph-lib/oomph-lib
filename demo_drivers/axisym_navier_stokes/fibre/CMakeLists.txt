# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered fibre subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(fibre C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME fibre
  SOURCES fibre.cc
  LIBRARIES oomph::axisym_navier_stokes oomph::fluid_interface)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME axisym_navier_stokes.fibre
  TARGET_DEPENDENCIES fibre
  EXTRA_REQUIRES
  LABELS axisym_navier_stokes fibre)

message(VERBOSE "Leaving fibre subdirectory")
# ------------------------------------------------------------------------------
