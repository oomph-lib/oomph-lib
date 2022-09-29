# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered double_diffusive_convection subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(double_diffusive_convection C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME dd_convection
  SOURCES dd_convection.cc db_navier_st_elements.h
  LIBRARIES oomph::advection_diffusion_reaction oomph::navier_stokes
  CXX_STANDARD 14)

oomph_add_executable(
  NAME multimesh_dd_convection
  SOURCES multimesh_dd_convection.cc
  LIBRARIES oomph::advection_diffusion oomph::navier_stokes
  CXX_STANDARD 14)

oomph_add_executable(
  NAME multimesh_ref_dd_convection
  SOURCES multimesh_ref_dd_convection.cc db_nst_external_elements.h
  LIBRARIES oomph::advection_diffusion oomph::navier_stokes
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME multi_physics.double_diffusive_convection
  TARGET_DEPENDENCIES dd_convection multimesh_dd_convection
                      multimesh_ref_dd_convection
  EXTRA_REQUIRES
  LABELS "multi_physics;double_diffusive_convection")

message(VERBOSE "Leaving double_diffusive_convection subdirectory")
# ------------------------------------------------------------------------------
