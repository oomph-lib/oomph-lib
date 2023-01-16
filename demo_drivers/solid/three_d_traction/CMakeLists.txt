# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered three_d_traction subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(three_d_traction C CXX Fortran)
if(NOT OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME three_d_traction
  SOURCES three_d_traction.cc
  LIBRARIES oomph::solid oomph::constitutive oomph::meshes oomph::generic)

oomph_add_executable(
  NAME three_d_traction_non_adapt
  SOURCES three_d_traction.cc
  LIBRARIES oomph::solid oomph::constitutive oomph::meshes oomph::generic
  CXX_DEFINITIONS NOREFINE)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME solid.three_d_traction
  TARGET_DEPENDENCIES three_d_traction three_d_traction_non_adapt
  VALIDATE_SH_ARGS ${OOMPH_ROOT_DIR}
  EXTRA_REQUIRES
  LABELS solid three_d_traction)

message(VERBOSE "Leaving three_d_traction subdirectory")
# ------------------------------------------------------------------------------
