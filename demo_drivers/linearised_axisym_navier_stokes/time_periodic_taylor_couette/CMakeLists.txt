# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered time_periodic_taylor_couette subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(time_periodic_taylor_couette C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME time_periodic_taylor_couette
  SOURCES ../self_starting_BDF2_timestepper.h time_periodic_taylor_couette.cc
  LIBRARIES oomph::axisym_navier_stokes oomph::linearised_axisym_navier_stokes
            oomph::navier_stokes oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME linearised_axisym_navier_stokes.time_periodic_taylor_couette
  TARGET_DEPENDENCIES time_periodic_taylor_couette
  EXTRA_REQUIRES
  LABELS linearised_axisym_navier_stokes time_periodic_taylor_couette)

message(VERBOSE "Leaving time_periodic_taylor_couette subdirectory")
# ------------------------------------------------------------------------------
