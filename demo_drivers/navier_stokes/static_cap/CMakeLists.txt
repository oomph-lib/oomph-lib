# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered static_cap subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(static_cap C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME static_single_layer
  SOURCES static_single_layer.cc
  LIBRARIES oomph::navier_stokes oomph::fluid_interface oomph::solid
            oomph::constitutive)

oomph_add_executable(
  NAME static_two_layer
  SOURCES static_two_layer.cc
  LIBRARIES oomph::navier_stokes oomph::fluid_interface oomph::solid
            oomph::constitutive)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME navier_stokes.static_cap
  TARGET_DEPENDENCIES static_single_layer static_two_layer
  EXTRA_REQUIRES
  LABELS navier_stokes static_cap)

message(VERBOSE "Leaving static_cap subdirectory")
# ------------------------------------------------------------------------------
