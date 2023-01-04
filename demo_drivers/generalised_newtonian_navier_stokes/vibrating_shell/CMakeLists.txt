# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered vibrating_shell subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(vibrating_shell C CXX Fortran)
if(NOT OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME vibrating_shell_non_newtonian
  SOURCES vibrating_shell_non_newtonian.cc overloaded_cartesian_element_body.h
  LIBRARIES oomph::generalised_newtonian_navier_stokes
            oomph::navier_stokes
            oomph::solid
            oomph::constitutive
            oomph::fluid_interface
            oomph::rigid_body
            oomph::meshes
            oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME generalised_newtonian_navier_stokes.vibrating_shell
  TARGET_DEPENDENCIES vibrating_shell_non_newtonian
  EXTRA_REQUIRES
  LABELS generalised_newtonian_navier_stokes vibrating_shell)

# This test needs more than the default allotted time to complete
set_tests_properties(generalised_newtonian_navier_stokes.vibrating_shell
                     PROPERTIES TIMEOUT 10000)

message(VERBOSE "Leaving vibrating_shell subdirectory")
# ------------------------------------------------------------------------------
