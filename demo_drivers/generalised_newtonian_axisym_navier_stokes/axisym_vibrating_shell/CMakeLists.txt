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
message(VERBOSE "Entered axisym_vibrating_shell subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
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
  DEPENDS_ON axisym_vibrating_shell_non_newtonian
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

# This test needs more than the default allotted time to complete
set_tests_properties(
  generalised_newtonian_axisym_navier_stokes.axisym_vibrating_shell
  PROPERTIES TIMEOUT 3600)

message(VERBOSE "Leaving axisym_vibrating_shell subdirectory")
# ------------------------------------------------------------------------------
