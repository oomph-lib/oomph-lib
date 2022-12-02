# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered counter_rotating_disks subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(counter_rotating_disks C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME counter_rotating_disks
  SOURCES ../self_starting_BDF2_timestepper.h counter_rotating_disks.cc
  LIBRARIES oomph::linearised_axisym_navier_stokes oomph::axisym_navier_stokes
            oomph::navier_stokes oomph::meshes oomph::generic)

oomph_add_executable(
  NAME counter_rotating_disks_ref
  SOURCES ../self_starting_BDF2_timestepper.h counter_rotating_disks_ref.cc
  LIBRARIES oomph::linearised_axisym_navier_stokes oomph::axisym_navier_stokes
            oomph::navier_stokes oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME linearised_axisym_navier_stokes.counter_rotating_disks
  TARGET_DEPENDENCIES counter_rotating_disks counter_rotating_disks_ref
  EXTRA_REQUIRES
  LABELS linearised_axisym_navier_stokes counter_rotating_disks)

# This test needs more than the default allotted time to complete
set_tests_properties(linearised_axisym_navier_stokes.counter_rotating_disks
                     PROPERTIES TIMEOUT 10000)

message(VERBOSE "Leaving counter_rotating_disks subdirectory")
# ------------------------------------------------------------------------------
