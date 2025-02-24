# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered complex_matrices_test subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
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

# TODO: Check this actually works and it doesn't just become an extra label
set(EXTRA_FLAGS)
if(OOMPH_HAS_MPI)
  set(EXTRA_FLAGS)
endif()

oomph_add_test(
  TEST_NAME self_test.generic.complex_matrices_test
  DEPENDS_ON dense_complex_matrix_test cr_complex_matrix_test
             cc_complex_matrix_test
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR} ${OOMPH_MPI_RUN_COMMAND}
  TEST_FILES validate.sh validata
  LABELS complex matrix ${EXTRA_FLAGS})

message(VERBOSE "Leaving complex_matrices_test subdirectory")
# ------------------------------------------------------------------------------
