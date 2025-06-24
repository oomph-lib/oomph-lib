# ------------------------------------------------------------------------------
# OOMPH auto-inserted: warn once if user configures inside an existing build.
if(EXISTS "${CMAKE_CURRENT_BINARY_DIR}/CMakeCache.txt"
   AND CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  message(
    WARNING
      "You're re-configuring an existing build tree. For a clean configuration, "
      "delete the whole build folder and start with a fresh directory.")
endif()
# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered counter_rotating_disks subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(counter_rotating_disks C CXX Fortran)

if(NOT oomphlib_FOUND)
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
  DEPENDS_ON counter_rotating_disks counter_rotating_disks_ref
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

# This test needs more than the default allotted time to complete
set_tests_properties(linearised_axisym_navier_stokes.counter_rotating_disks
                     PROPERTIES TIMEOUT 3600)

message(VERBOSE "Leaving counter_rotating_disks subdirectory")
# ------------------------------------------------------------------------------
