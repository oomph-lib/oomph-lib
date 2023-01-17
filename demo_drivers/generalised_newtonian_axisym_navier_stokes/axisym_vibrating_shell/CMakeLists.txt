# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered axisym_vibrating_shell subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(axisym_vibrating_shell C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME axisym_vibrating_shell_non_newtonian
  SOURCES axisym_vibrating_shell_non_newtonian.cc overloaded_element_body.h
  LIBRARIES oomph::generalised_newtonian_axisym_navier_stokes
            oomph::axisym_navier_stokes
            oomph::solid
            oomph::constitutive
            oomph::fluid_interface
            oomph::rigid_body
            oomph::meshes
            oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME generalised_newtonian_axisym_navier_stokes.axisym_vibrating_shell
  TARGET_DEPENDENCIES axisym_vibrating_shell_non_newtonian
  VALIDATE_SH_ARGS ${OOMPH_ROOT_DIR}
  EXTRA_REQUIRES
  LABELS generalised_newtonian_axisym_navier_stokes axisym_vibrating_shell)

# This test needs more than the default allotted time to complete
set_tests_properties(
  generalised_newtonian_axisym_navier_stokes.axisym_vibrating_shell
  PROPERTIES TIMEOUT 10000)

message(VERBOSE "Leaving axisym_vibrating_shell subdirectory")
# ------------------------------------------------------------------------------
