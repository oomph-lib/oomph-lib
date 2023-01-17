# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered point_source subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(point_source C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME adaptive_helmholtz_point_source
  SOURCES helmholtz_point_source.cc
  LIBRARIES oomph::helmholtz oomph::meshes oomph::generic
  CXX_DEFINITIONS ADAPTIVE)

oomph_add_executable(
  NAME helmholtz_point_source
  SOURCES helmholtz_point_source.cc
  LIBRARIES oomph::helmholtz oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME helmholtz.point_source
  TARGET_DEPENDENCIES adaptive_helmholtz_point_source helmholtz_point_source
  VALIDATE_SH_ARGS ${OOMPH_ROOT_DIR}
  EXTRA_REQUIRES
  LABELS helmholtz point_source)

message(VERBOSE "Leaving point_source subdirectory")
# ------------------------------------------------------------------------------
