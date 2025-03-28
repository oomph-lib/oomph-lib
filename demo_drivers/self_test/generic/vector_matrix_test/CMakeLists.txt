# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered vector_matrix_test subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(vector_matrix_test C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME vector_matrix_test
  SOURCES vector_matrix_test.cc
  LIBRARIES oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME self_test.generic.vector_matrix_test
  DEPENDS_ON vector_matrix_test
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving vector_matrix_test subdirectory")
# ------------------------------------------------------------------------------
