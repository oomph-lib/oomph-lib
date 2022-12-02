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
  SOURCES counter_rotating_disks.cc linearised_axisym_navier_stokes_elements.h
          multi_domain_linearised_axisym_navier_stokes_elements.h
  LIBRARIES oomph::axisym_navier_stokes oomph::navier_stokes oomph::meshes
            oomph::generic)

oomph_add_executable(
  NAME counter_rotating_disks_ref
  SOURCES counter_rotating_disks_ref.cc
          refineable_linearised_axisym_navier_stokes_elements.h
          linearised_axisym_navier_stokes_elements.h
          multi_domain_linearised_axisym_navier_stokes_elements.h
  LIBRARIES oomph::axisym_navier_stokes oomph::navier_stokes oomph::meshes
            oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME axisym_navier_stokes.counter_rotating_disks
  TARGET_DEPENDENCIES counter_rotating_disks counter_rotating_disks_ref
  EXTRA_REQUIRES
  LABELS axisym_navier_stokes counter_rotating_disks)

message(VERBOSE "Leaving counter_rotating_disks subdirectory")
# ------------------------------------------------------------------------------
