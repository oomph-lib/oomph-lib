# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered unstructured_two_d_circle subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(unstructured_two_d_circle C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME adaptive_unstructured_two_d_circle
  SOURCES unstructured_two_d_circle.cc
  LIBRARIES oomph::darcy oomph::meshes oomph::generic
  CXX_DEFINITIONS ADAPTIVE)

oomph_add_executable(
  NAME unstructured_two_d_circle
  SOURCES unstructured_two_d_circle.cc
  LIBRARIES oomph::darcy oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME darcy.unstructured_two_d_circle
  TARGET_DEPENDENCIES adaptive_unstructured_two_d_circle
                      unstructured_two_d_circle
  VALIDATE_SH_ARGS ${OOMPH_ROOT_DIR}
  EXTRA_REQUIRES
  LABELS darcy unstructured_two_d_circle)

message(VERBOSE "Leaving unstructured_two_d_circle subdirectory")
# ------------------------------------------------------------------------------
