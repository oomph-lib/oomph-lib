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
message(VERBOSE "Entered airy_cantilever subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(airy_cantilever C CXX Fortran)

if(NOT oomphlib_FOUND)
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
  TEST_NAME mpi.distribution.airy_cantilever
  DEPENDS_ON airy_cantilever2_adapt
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR} ${OOMPH_MPI_RUN_COMMAND}
  TEST_FILES validate.sh
             validata
             airy_cantilever_0_partition.dat
             airy_cantilever_1_partition.dat
             airy_cantilever_2_partition.dat
             airy_cantilever_3_partition.dat
             airy_cantilever_4_partition.dat
             airy_cantilever_5_partition.dat
             airy_cantilever_6_partition.dat
             airy_cantilever_7_partition.dat
             airy_cantilever_8_partition.dat
             airy_cantilever_9_partition.dat)

# Don't run this test alongside any other tests
set_tests_properties(mpi.distribution.airy_cantilever
                     PROPERTIES PROCESSORS ${OOMPH_MPI_NUM_PROC})

message(VERBOSE "Leaving airy_cantilever subdirectory")
# ------------------------------------------------------------------------------
