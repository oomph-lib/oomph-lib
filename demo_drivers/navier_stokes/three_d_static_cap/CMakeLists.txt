# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered three_d_static_cap subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(three_d_static_cap C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME 3d_static_cap
  SOURCES 3d_static_cap.cc
  LIBRARIES oomph::navier_stokes oomph::fluid_interface
  CXX_STANDARD 17)

oomph_add_executable(
  NAME 3d_static_cap_elastic
  SOURCES 3d_static_cap_elastic.cc
  LIBRARIES oomph::navier_stokes oomph::fluid_interface oomph::solid
            oomph::constitutive
  CXX_STANDARD 17)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME navier_stokes.three_d_static_cap
  TARGET_DEPENDENCIES 3d_static_cap 3d_static_cap_elastic
  EXTRA_REQUIRES
  LABELS navier_stokes three_d_static_cap)

message(VERBOSE "Leaving three_d_static_cap subdirectory")
# ------------------------------------------------------------------------------
