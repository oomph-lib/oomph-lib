# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered eigen_solver subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(eigen_solver C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME eigen_solver_test
  SOURCES eigen_solver_test.cc
  LIBRARIES oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME mpi.eigen_solver
  DEPENDS_ON eigen_solver
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata random_test_matrix.dat
  LABELS mpi matrix)

message(VERBOSE "Leaving eigen_solver subdirectory")
# ------------------------------------------------------------------------------
