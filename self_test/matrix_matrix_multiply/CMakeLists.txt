# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered matrix_matrix_multiply subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(matrix_matrix_multiply C CXX Fortran)
if(NOT OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME matrix_multiply_test
  SOURCES matrix_multiply_test.cc validate.sh
  LIBRARIES oomph::navier_stokes oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME matrix_matrix_multiply
  TARGET_DEPENDENCIES matrix_multiply_test
  EXTRA_REQUIRES
  LABELS matrix)

message(VERBOSE "Leaving matrix_matrix_multiply subdirectory")
# ------------------------------------------------------------------------------
