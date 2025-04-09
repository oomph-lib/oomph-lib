# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered problem_test subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(problem_test C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME solve_eigenproblem_test
  SOURCES solve_eigenproblem_test.cc
  LIBRARIES oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME self_test.generic.problem_test
  DEPENDS_ON solve_eigenproblem_test
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata random_test_matrix.dat)

message(VERBOSE "Leaving problem_test subdirectory")
# ------------------------------------------------------------------------------
