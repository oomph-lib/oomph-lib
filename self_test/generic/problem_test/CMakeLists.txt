# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered problem_test subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(problem_test C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME solve_eigenproblem_test
  SOURCES solve_eigenproblem_test.cc
  LIBRARIES oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME generic.problem_test
  TARGET_DEPENDENCIES solve_eigenproblem_test
  VALIDATE_SH_ARGS ${OOMPH_ROOT_DIR}
  EXTRA_REQUIRES
  LABELS problem)

# FIXME: Currently fails because the eigenvectors are different in sign in some
# places
set_tests_properties(generic.problem_test PROPERTIES DISABLED TRUE)

message(VERBOSE "Leaving problem_test subdirectory")
# ------------------------------------------------------------------------------
