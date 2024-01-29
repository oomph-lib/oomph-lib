# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered eigen_solver_test subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(eigen_solver_test C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME eigen_solver_test
  SOURCES eigen_solver_test.cc
  LIBRARIES oomph::meshes oomph::generic)

oomph_add_executable(
  NAME lapack_qz_find_eigenvalues_test
  SOURCES lapack_qz_find_eigenvalues_test.cc
  LIBRARIES oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME self_test.generic.eigen_solver_test
  DEPENDS_ON eigen_solver_test lapack_qz_find_eigenvalues_test
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata random_test_matrix.dat
  LABELS eigen_solver lapack)

# FIXME: Solution is completely off...
set_tests_properties(self_test.generic.eigen_solver_test PROPERTIES DISABLED
                                                                    TRUE)

message(VERBOSE "Leaving eigen_solver_test subdirectory")
# ------------------------------------------------------------------------------
