# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered line_visualiser subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(line_visualiser C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME adaptive_driven_cavity
  SOURCES adaptive_driven_cavity.cc
  LIBRARIES oomph::navier_stokes oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME self_test.mpi.line_visualiser
  DEPENDS_ON adaptive_driven_cavity
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR} ${OOMPH_MPI_RUN_COMMAND}
  TEST_FILES validate.sh validata adaptive_cavity_1_partition.dat
             adaptive_cavity_2_partition.dat
  LABELS mpi adapt driven_cavity)

# Don't run this test alongside any other tests
set_tests_properties(self_test.mpi.line_visualiser
                     PROPERTIES PROCESSORS ${OOMPH_MPI_NUM_PROC})

message(VERBOSE "Leaving line_visualiser subdirectory")
# ------------------------------------------------------------------------------
