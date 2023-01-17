# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered C_style_output subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
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
  TARGET_DEPENDENCIES c_style_output
  VALIDATE_SH_ARGS ${OOMPH_ROOT_DIR}
  EXTRA_REQUIRES
  LABELS optimisation C_style_output
  SILENCE_MISSING_VALIDATA_WARNING)

message(VERBOSE "Leaving C_style_output subdirectory")
# ------------------------------------------------------------------------------
