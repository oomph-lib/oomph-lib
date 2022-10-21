# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered linear_solvers subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(linear_solvers C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME adv_diff_iterative_linear_solver_tester
  SOURCES adv_diff_iterative_linear_solver_tester.cc
  LIBRARIES oomph::advection_diffusion
  CXX_STANDARD 14)

oomph_add_executable(
  NAME driven_cavity
  SOURCES driven_cavity.cc
  LIBRARIES oomph::navier_stokes
  CXX_STANDARD 14)

oomph_add_executable(
  NAME direct_solver_test
  SOURCES direct_solver_test.cc
  LIBRARIES oomph::poisson
  CXX_STANDARD 14)

oomph_add_executable(
  NAME two_d_poisson_dense_lu
  SOURCES two_d_poisson_dense_lu.cc
  LIBRARIES oomph::poisson
  CXX_STANDARD 14)

oomph_add_executable(
  NAME two_d_linear_elasticity_with_simple_block_diagonal_preconditioner
  SOURCES two_d_linear_elasticity_with_simple_block_diagonal_preconditioner.cc
          simple_block_preconditioners.h
  LIBRARIES oomph::linear_elasticity
  CXX_STANDARD 14)

oomph_add_executable(
  NAME two_d_multi_poisson
  SOURCES two_d_multi_poisson.cc multi_poisson_block_preconditioners.h
          multi_poisson_elements.h
  LIBRARIES oomph::generic
  CXX_STANDARD 14)

oomph_add_executable(
  NAME two_d_poisson_cg
  SOURCES two_d_poisson_cg.cc
  LIBRARIES oomph::poisson
  CXX_STANDARD 14)

# Don't add these drivers if we have MPI as Hypre and Trilinos can't handle
# serial execution of MPI-enabled code
if(NOT OOMPH_HAS_MPI)
  if(OOMPH_HAS_HYPRE)
    oomph_add_executable(
      NAME HypreSolver_test
      SOURCES HypreSolver_test.cc
      LIBRARIES oomph::poisson
      CXX_STANDARD 14)
  endif()

  if(OOMPH_HAS_TRILINOS)
    oomph_add_executable(
      NAME TrilinosSolver_test
      SOURCES TrilinosSolver_test.cc
      LIBRARIES oomph::poisson Trilinos::Trilinos
      CXX_STANDARD 14)
  endif()
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

oomph_add_test(
  TEST_NAME linear_solvers.linear_solvers
  TARGET_DEPENDENCIES ${TEST_TARGET_DEPS}
  EXTRA_REQUIRES
  LABELS "linear_solvers")

message(VERBOSE "Leaving linear_solvers subdirectory")
# ------------------------------------------------------------------------------
