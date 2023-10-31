# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered falling_jet subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(falling_jet C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME falling_jet
  SOURCES falling_jet.cc
  LIBRARIES oomph::navier_stokes oomph::fluid_interface oomph::constitutive
            oomph::solid oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME navier_stokes.falling_jet
  DEPENDS_ON falling_jet
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata
  LABELS navier_stokes falling_jet)

message(VERBOSE "Leaving falling_jet subdirectory")
# ------------------------------------------------------------------------------
