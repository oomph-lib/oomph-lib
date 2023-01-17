# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered vector_matrix_test subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(vector_matrix_test C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME vector_matrix_test
  SOURCES vector_matrix_test.cc
  LIBRARIES oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME generic.vector_matrix_test
  TARGET_DEPENDENCIES vector_matrix_test
  VALIDATE_SH_ARGS ${OOMPH_ROOT_DIR}
  EXTRA_REQUIRES
  LABELS vector matrix)

message(VERBOSE "Leaving vector_matrix_test subdirectory")
# ------------------------------------------------------------------------------
