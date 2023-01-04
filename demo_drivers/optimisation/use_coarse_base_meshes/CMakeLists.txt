# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered use_coarse_base_meshes subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(use_coarse_base_meshes C CXX Fortran)
if(NOT OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME two_d_poisson_adapt
  SOURCES two_d_poisson_adapt.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME optimisation.use_coarse_base_meshes
  TARGET_DEPENDENCIES two_d_poisson_adapt
  EXTRA_REQUIRES
  LABELS optimisation use_coarse_base_meshes)

message(VERBOSE "Leaving use_coarse_base_meshes subdirectory")
# ------------------------------------------------------------------------------
