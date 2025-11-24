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
message(VERBOSE "Entered linear_solvers subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(linear_solvers C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../install")
endif()

include(CTest)

oomph_add_executable(
  NAME adv_diff_iterative_linear_solver_tester
  SOURCES adv_diff_iterative_linear_solver_tester.cc
  LIBRARIES oomph::advection_diffusion oomph::meshes oomph::generic)

oomph_add_executable(
  NAME driven_cavity
  SOURCES driven_cavity.cc
  LIBRARIES oomph::navier_stokes oomph::meshes oomph::generic)

oomph_add_executable(
  NAME direct_solver_test
  SOURCES direct_solver_test.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

oomph_add_executable(
  NAME two_d_poisson_dense_lu
  SOURCES two_d_poisson_dense_lu.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

oomph_add_executable(
  NAME two_d_linear_elasticity_with_simple_block_diagonal_preconditioner
  SOURCES two_d_linear_elasticity_with_simple_block_diagonal_preconditioner.cc
          simple_block_preconditioners.h
  LIBRARIES oomph::linear_elasticity oomph::meshes oomph::generic)

oomph_add_executable(
  NAME two_d_multi_poisson
  SOURCES two_d_multi_poisson.cc multi_poisson_block_preconditioners.h
          multi_poisson_elements.h
  LIBRARIES oomph::meshes oomph::generic)

oomph_add_executable(
  NAME two_d_poisson_cg
  SOURCES two_d_poisson_cg.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

# Don't add these drivers if we have MPI as Hypre and Trilinos can't handle
# serial execution of MPI-enabled code
if(NOT OOMPH_HAS_MPI)
  if(OOMPH_HAS_HYPRE)
    oomph_add_executable(
      NAME HypreSolver_test
      SOURCES HypreSolver_test.cc
      LIBRARIES oomph::poisson oomph::meshes oomph::generic)
  endif()

  if(OOMPH_HAS_TRILINOS)
    oomph_add_executable(
      NAME TrilinosSolver_test
      SOURCES TrilinosSolver_test.cc
      LIBRARIES oomph::poisson oomph::meshes oomph::generic)
  endif()
endif()

if(OOMPH_HAS_MUMPS)
  oomph_add_executable(
    NAME serial_mumps_test
    SOURCES serial_mumps_test.cc
    LIBRARIES oomph::poisson oomph::meshes oomph::generic)
endif()

# ---------------------------------- TESTING ----------------------------------

set(TEST_TARGET_DEPS
    adv_diff_iterative_linear_solver_tester
    driven_cavity
    direct_solver_test
    two_d_poisson_dense_lu
    two_d_linear_elasticity_with_simple_block_diagonal_preconditioner
    two_d_multi_poisson
    two_d_poisson_cg)

if(NOT OOMPH_HAS_MPI)
  if(OOMPH_HAS_HYPRE)
    list(APPEND TEST_TARGET_DEPS HypreSolver_test)
  endif()

  if(OOMPH_HAS_TRILINOS)
    list(APPEND TEST_TARGET_DEPS TrilinosSolver_test)
  endif()
endif()

if(OOMPH_HAS_MUMPS)
  list(APPEND TEST_TARGET_DEPS serial_mumps_test)
endif()

oomph_add_test(
  TEST_NAME linear_solvers.linear_solvers
  DEPENDS_ON ${TEST_TARGET_DEPS}
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving linear_solvers subdirectory")
# ------------------------------------------------------------------------------
