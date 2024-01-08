# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered boussinesq_convection subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(boussinesq_convection C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME multi_domain_boussinesq_convection
  SOURCES multi_domain_boussinesq_convection.cc
  LIBRARIES oomph::advection_diffusion oomph::navier_stokes oomph::meshes
            oomph::generic)

oomph_add_executable(
  NAME multi_domain_boussinesq_convection_fd
  SOURCES multi_domain_boussinesq_convection.cc
  LIBRARIES oomph::advection_diffusion oomph::navier_stokes oomph::meshes
            oomph::generic
  CXX_DEFINITIONS USE_FD_JACOBIAN_NST_IN_MULTI_DOMAIN_BOUSSINESQ)

oomph_add_executable(
  NAME multi_domain_ref_b_convection
  SOURCES multi_domain_ref_b_convection.cc
  LIBRARIES oomph::advection_diffusion oomph::navier_stokes oomph::meshes
            oomph::generic)

oomph_add_executable(
  NAME refineable_b_convection
  SOURCES refineable_b_convection.cc
  LIBRARIES oomph::advection_diffusion oomph::navier_stokes oomph::meshes
            oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME mpi.boussinesq_convection
  DEPENDS_ON multi_domain_boussinesq_convection
             multi_domain_boussinesq_convection_fd multi_domain_ref_b_convection
             refineable_b_convection
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR} ${OOMPH_MPI_RUN_COMMAND}
  TEST_FILES validate.sh
             validata
             multi_domain_boussinesq_partition.dat
             multi_domain_ref_b_partition.dat
             refineable_b_convection_partition.dat
             multi_domain_boussinesq_partition_2.dat
             prune_and_balance.mcr
  LABELS mpi boussinesq_convection)

# Don't run this test alongside any other tests
set_tests_properties(mpi.boussinesq_convection PROPERTIES PROCESSORS
                                                          ${OOMPH_MPI_NUM_PROC})
if(APPLE)
  set_tests_properties(mpi.boussinesq_convection PROPERTIES DISABLED FALSE)
endif()

message(VERBOSE "Leaving boussinesq_convection subdirectory")
# ------------------------------------------------------------------------------
