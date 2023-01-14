# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered eigen_solver subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(eigen_solver C CXX Fortran)
if(NOT OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY)
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
  TARGET_DEPENDENCIES eigen_solver
  VALIDATE_SH_ARGS ${OOMPH_ROOT_DIR}
  EXTRA_REQUIRES random_test_matrix.dat
  LABELS mpi matrix)

message(VERBOSE "Leaving eigen_solver subdirectory")
# ------------------------------------------------------------------------------
