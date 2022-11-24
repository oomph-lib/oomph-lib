# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered 3d_boussinesq_convection subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(3d_boussinesq_convection C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME 3d_multimesh_ref_b_convect
  SOURCES 3d_multimesh_ref_b_convect.cc my_boussinesq_elements.h
  LIBRARIES oomph::advection_diffusion oomph::navier_stokes)

oomph_add_executable(
  NAME 3d_ref_b_convect
  SOURCES 3d_ref_b_convect.cc
  LIBRARIES oomph::advection_diffusion oomph::navier_stokes)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME multi_physics.3d_boussinesq_convection
  TARGET_DEPENDENCIES 3d_multimesh_ref_b_convect 3d_ref_b_convect
  EXTRA_REQUIRES
  LABELS "multi_physics;3d_boussinesq_convection")

message(VERBOSE "Leaving 3d_boussinesq_convection subdirectory")
# ------------------------------------------------------------------------------
