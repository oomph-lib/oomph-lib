# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered boussinesq_convection subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(boussinesq_convection C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME boussinesq_convection
  SOURCES boussinesq_convection.cc
  LIBRARIES oomph::advection_diffusion oomph::navier_stokes
  CXX_DEFINITIONS USE_OFF_DIAGONAL_FD_JACOBIAN_FOR_BUOYANT_Q_ELEMENT
  CXX_STANDARD 17)

oomph_add_executable(
  NAME boussinesq_convection_analytic
  SOURCES boussinesq_convection.cc
  LIBRARIES oomph::advection_diffusion oomph::navier_stokes
  CXX_STANDARD 17)

oomph_add_executable(
  NAME boussinesq_convection_fd
  SOURCES boussinesq_convection.cc
  LIBRARIES oomph::advection_diffusion oomph::navier_stokes
  CXX_DEFINITIONS USE_FD_JACOBIAN_FOR_BUOYANT_Q_ELEMENT
  CXX_STANDARD 17)

oomph_add_executable(
  NAME multi_domain_boussinesq_convection
  SOURCES multi_domain_boussinesq_convection.cc
  LIBRARIES oomph::advection_diffusion oomph::navier_stokes
  CXX_STANDARD 17)

oomph_add_executable(
  NAME multi_domain_ref_b_convection
  SOURCES multi_domain_ref_b_convection.cc
  LIBRARIES oomph::advection_diffusion oomph::navier_stokes
  CXX_STANDARD 17)

oomph_add_executable(
  NAME multi_domain_ref_b_convection_fd_for_external_data
  SOURCES multi_domain_ref_b_convection.cc
  LIBRARIES oomph::advection_diffusion oomph::navier_stokes
  CXX_DEFINITIONS
    USE_FD_FOR_DERIVATIVES_WRT_EXTERNAL_DATA_IN_MULTI_DOMAIN_BOUSSINESQ
  CXX_STANDARD 17)

oomph_add_executable(
  NAME refineable_b_convection
  SOURCES refineable_b_convection.cc
  LIBRARIES oomph::advection_diffusion oomph::navier_stokes
  CXX_STANDARD 17)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME multi_physics.boussinesq_convection
  TARGET_DEPENDENCIES
    boussinesq_convection
    boussinesq_convection_analytic
    boussinesq_convection_fd
    multi_domain_boussinesq_convection
    multi_domain_ref_b_convection
    multi_domain_ref_b_convection_fd_for_external_data
    refineable_b_convection
  EXTRA_REQUIRES
  LABELS multi_physics boussinesq_convection)

message(VERBOSE "Leaving boussinesq_convection subdirectory")
# ------------------------------------------------------------------------------
