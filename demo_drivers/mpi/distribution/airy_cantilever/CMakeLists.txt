# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered airy_cantilever subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(airy_cantilever C CXX Fortran)
if(NOT OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME airy_cantilever2_adapt
  SOURCES airy_cantilever2.cc
  LIBRARIES oomph::solid oomph::constitutive oomph::meshes oomph::generic
  CXX_DEFINITIONS REFINE)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME mpi.airy_cantilever
  TARGET_DEPENDENCIES airy_cantilever2_adapt
  VALIDATE_SH_ARGS ${OOMPH_ROOT_DIR} ${OOMPH_MPI_RUN_COMMAND}
  EXTRA_REQUIRES
    airy_cantilever_0_partition.dat
    airy_cantilever_1_partition.dat
    airy_cantilever_2_partition.dat
    airy_cantilever_3_partition.dat
    airy_cantilever_4_partition.dat
    airy_cantilever_5_partition.dat
    airy_cantilever_6_partition.dat
    airy_cantilever_7_partition.dat
    airy_cantilever_8_partition.dat
    airy_cantilever_9_partition.dat
  LABELS mpi airy_cantilever)

# Don't run this test alongside any other tests
set_tests_properties(mpi.airy_cantilever PROPERTIES PROCESSORS
                                                    ${OOMPH_MPI_NUM_PROC})

message(VERBOSE "Leaving airy_cantilever subdirectory")
# ------------------------------------------------------------------------------
