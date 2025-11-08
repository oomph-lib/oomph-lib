# ------------------------------------------------------------------------------
# OOMPH auto-inserted: warn once if user configures inside an existing build.
if(EXISTS "${CMAKE_CURRENT_BINARY_DIR}/CMakeCache.txt"
   AND CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  message(
    WARNING
      "You're re-configuring an existing build tree. For a clean configuration, "
      "delete the whole build folder and start with a fresh directory.")
endif()
# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered complex_matrices_test subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(complex_matrices_test C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../../install")
endif()

include(CTest)

oomph_add_executable(
  NAME dense_complex_matrix_test
  SOURCES dense_complex_matrix_test.cc
  LIBRARIES oomph::generic)

oomph_add_executable(
  NAME cr_complex_matrix_test
  SOURCES cr_complex_matrix_test.cc
  LIBRARIES oomph::generic)

oomph_add_executable(
  NAME cc_complex_matrix_test
  SOURCES cc_complex_matrix_test.cc
  LIBRARIES oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME self_test.generic.complex_matrices_test
  DEPENDS_ON dense_complex_matrix_test cr_complex_matrix_test
             cc_complex_matrix_test
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR} ${OOMPH_MPI_RUN_COMMAND}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving complex_matrices_test subdirectory")
# ------------------------------------------------------------------------------
