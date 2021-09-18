# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered counter_rotating_disks subdirectory")

cmake_minimum_required(VERSION 3.16..3.20 FATAL_ERROR)
project(counter_rotating_disks C CXX Fortran)
if(NOT BUILDING_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME counter_rotating_disks
  SOURCES counter_rotating_disks.cc linearised_axisym_navier_stokes_elements.h
          multi_domain_linearised_axisym_navier_stokes_elements.h
  LIBRARIES oomph::axisym_navier_stokes oomph::navier_stokes
  CXX_STANDARD 14)

oomph_add_executable(
  NAME counter_rotating_disks_ref
  SOURCES counter_rotating_disks_ref.cc
          refineable_linearised_axisym_navier_stokes_elements.h
          linearised_axisym_navier_stokes_elements.h
          multi_domain_linearised_axisym_navier_stokes_elements.h
  LIBRARIES oomph::axisym_navier_stokes oomph::navier_stokes
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME axisym_navier_stokes.counter_rotating_disks
  TARGET_DEPENDENCIES counter_rotating_disks counter_rotating_disks_ref
  EXTRA_REQUIRES
  LABELS "axisym_navier_stokes;counter_rotating_disks")

message(VERBOSE "Leaving counter_rotating_disks subdirectory")
# ------------------------------------------------------------------------------
