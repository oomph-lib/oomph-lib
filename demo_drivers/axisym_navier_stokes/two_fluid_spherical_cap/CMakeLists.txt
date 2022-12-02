# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered two_fluid_spherical_cap subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(two_fluid_spherical_cap C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME two_fluid_spherical_cap
  SOURCES two_fluid_spherical_cap.cc
  LIBRARIES oomph::axisym_navier_stokes oomph::fluid_interface oomph::meshes
            oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME axisym_navier_stokes.two_fluid_spherical_cap
  TARGET_DEPENDENCIES two_fluid_spherical_cap
  EXTRA_REQUIRES
  LABELS axisym_navier_stokes two_fluid_spherical_cap)

message(VERBOSE "Leaving two_fluid_spherical_cap subdirectory")
# ------------------------------------------------------------------------------
