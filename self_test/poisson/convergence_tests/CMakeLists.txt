# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered convergence_tests subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(convergence_tests C CXX Fortran)
if(NOT OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME q_convergence_3d
  SOURCES q_convergence_3d.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

oomph_add_executable(
  NAME q_convergence_2d
  SOURCES q_convergence_2d.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

oomph_add_executable(
  NAME t_convergence_2d
  SOURCES t_convergence_2d.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

oomph_add_executable(
  NAME t_convergence_3d
  SOURCES t_convergence_3d.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME poisson.convergence_tests
  TARGET_DEPENDENCIES q_convergence_3d q_convergence_2d t_convergence_2d
                      t_convergence_3d
  EXTRA_REQUIRES
  LABELS convergence poisson)

message(VERBOSE "Leaving convergence_tests subdirectory")
# ------------------------------------------------------------------------------
