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
message(VERBOSE "Entered C_style_output subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(C_style_output C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()

include(CTest)

oomph_add_executable(
  NAME c_style_output
  SOURCES c_style_output.cc
  LIBRARIES oomph::unsteady_heat
            oomph::navier_stokes
            oomph::poisson
            oomph::beam
            oomph::solid
            oomph::constitutive
            oomph::meshes
            oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME optimisation.C_style_output
  DEPENDS_ON c_style_output
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh
  SILENCE_MISSING_VALIDATA_WARNING)

message(VERBOSE "Leaving C_style_output subdirectory")
# ------------------------------------------------------------------------------
