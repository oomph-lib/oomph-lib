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
message(VERBOSE "Entered boussinesq_convection subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(boussinesq_convection C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()

include(CTest)

oomph_add_executable(
  NAME boussinesq_convection
  SOURCES boussinesq_convection.cc
  LIBRARIES oomph::advection_diffusion oomph::navier_stokes oomph::meshes
            oomph::generic
  CXX_DEFINITIONS USE_OFF_DIAGONAL_FD_JACOBIAN_FOR_BUOYANT_Q_ELEMENT)

oomph_add_executable(
  NAME boussinesq_convection_analytic
  SOURCES boussinesq_convection.cc
  LIBRARIES oomph::advection_diffusion oomph::navier_stokes oomph::meshes
            oomph::generic)

oomph_add_executable(
  NAME boussinesq_convection_fd
  SOURCES boussinesq_convection.cc
  LIBRARIES oomph::advection_diffusion oomph::navier_stokes oomph::meshes
            oomph::generic
  CXX_DEFINITIONS USE_FD_JACOBIAN_FOR_BUOYANT_Q_ELEMENT)

oomph_add_executable(
  NAME multi_domain_boussinesq_convection
  SOURCES multi_domain_boussinesq_convection.cc
  LIBRARIES oomph::advection_diffusion oomph::navier_stokes oomph::meshes
            oomph::generic)

oomph_add_executable(
  NAME multi_domain_ref_b_convection
  SOURCES multi_domain_ref_b_convection.cc
  LIBRARIES oomph::advection_diffusion oomph::navier_stokes oomph::meshes
            oomph::generic)

oomph_add_executable(
  NAME multi_domain_ref_b_convection_fd_for_external_data
  SOURCES multi_domain_ref_b_convection.cc
  LIBRARIES oomph::advection_diffusion oomph::navier_stokes oomph::meshes
            oomph::generic
  CXX_DEFINITIONS
    USE_FD_FOR_DERIVATIVES_WRT_EXTERNAL_DATA_IN_MULTI_DOMAIN_BOUSSINESQ)

oomph_add_executable(
  NAME refineable_b_convection
  SOURCES refineable_b_convection.cc
  LIBRARIES oomph::advection_diffusion oomph::navier_stokes oomph::meshes
            oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME multi_physics.boussinesq_convection
  DEPENDS_ON boussinesq_convection
             boussinesq_convection_analytic
             boussinesq_convection_fd
             multi_domain_boussinesq_convection
             multi_domain_ref_b_convection
             multi_domain_ref_b_convection_fd_for_external_data
             refineable_b_convection
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving boussinesq_convection subdirectory")
# ------------------------------------------------------------------------------
