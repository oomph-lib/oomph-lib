# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered steady_rot subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(steady_rot C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME steady_rot
  SOURCES steady_rot.cc
  LIBRARIES oomph::spherical_navier_stokes oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME spherical_navier_stokes.steady_rot
  TARGET_DEPENDENCIES steady_rot
  VALIDATE_SH_ARGS ${OOMPH_ROOT_DIR}
  EXTRA_REQUIRES
  LABELS spherical_navier_stokes steady_rot)

message(VERBOSE "Leaving steady_rot subdirectory")
# ------------------------------------------------------------------------------
