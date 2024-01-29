# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered matrix_matrix_multiply subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(matrix_matrix_multiply C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME matrix_multiply_test
  SOURCES matrix_multiply_test.cc validate.sh
  LIBRARIES oomph::navier_stokes oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME self_test.matrix_matrix_multiply
  DEPENDS_ON matrix_multiply_test
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata
  LABELS matrix)

message(VERBOSE "Leaving matrix_matrix_multiply subdirectory")
# ------------------------------------------------------------------------------
