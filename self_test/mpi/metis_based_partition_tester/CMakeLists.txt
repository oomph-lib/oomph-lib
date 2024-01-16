# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered metis_based_partition_tester subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(metis_based_partition_tester C CXX Fortran)
if(NOT oomphlib_FOUND)
  # hierher the path to install needs to be updated on a case by case basis,
  # right?
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME adaptive_driven_cavity
  SOURCES adaptive_driven_cavity.cc
  LIBRARIES oomph::navier_stokes oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

# hierher: Puneet: I'm confused about the test name; I thought this was set
# automatically. If not, should we explain our conventions here? Can I do
# anything I want? Same with labels. Is there a master list?
oomph_add_test(
  TEST_NAME self_test.mpi.metis_based_partition_tester
  DEPENDS_ON adaptive_driven_cavity
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR} ${OOMPH_MPI_RUN_COMMAND}
  TEST_FILES validate.sh validata
  LABELS self_test mpi adaptive_driven_cavity)

# Don't run this test alongside any other tests hierher Puneet: how does ensure
# that it's not run alongside other tests? Does the first arg have to match the
# test name?
set_tests_properties(self_test.mpi.metis_based_partition_tester
                     PROPERTIES PROCESSORS ${OOMPH_MPI_NUM_PROC})

message(VERBOSE "Leaving metis_based_partition_tester subdirectory")
# ------------------------------------------------------------------------------
