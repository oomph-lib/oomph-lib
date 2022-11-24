# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered jeffery_orbit subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(jeffery_orbit C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME jeffery_orbit
  SOURCES jeffery_orbit.cc my_taylor_hood_elements.h
  LIBRARIES oomph::solid oomph::constitutive oomph::navier_stokes
            oomph::rigid_body)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME navier_stokes.jeffery_orbit
  TARGET_DEPENDENCIES jeffery_orbit
  EXTRA_REQUIRES
  LABELS navier_stokes jeffery_orbit)

message(VERBOSE "Leaving jeffery_orbit subdirectory")
# ------------------------------------------------------------------------------
