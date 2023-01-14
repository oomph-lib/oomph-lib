# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered complex_matrices_test subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(complex_matrices_test C CXX Fortran)
if(NOT OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
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
  TEST_NAME generic.complex_matrices_test
  TARGET_DEPENDENCIES dense_complex_matrix_test cr_complex_matrix_test
                      cc_complex_matrix_test
  VALIDATE_SH_ARGS ${OOMPH_ROOT_DIR} ${OOMPH_MPI_RUN_COMMAND}
  EXTRA_REQUIRES
  LABELS complex matrix ${EXTRA_FLAGS})

set_tests_properties(generic.complex_matrices_test PROPERTIES DISABLED TRUE)

message(VERBOSE "Leaving complex_matrices_test subdirectory")
# ------------------------------------------------------------------------------
