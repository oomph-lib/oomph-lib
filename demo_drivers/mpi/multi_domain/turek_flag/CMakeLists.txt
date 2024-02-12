# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered turek_flag subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(turek_flag C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME turek_flag
  SOURCES turek_flag.cc
  LIBRARIES oomph::solid oomph::constitutive oomph::navier_stokes oomph::meshes
            oomph::generic)

oomph_add_executable(
  NAME turek_flag_load_balance
  SOURCES turek_flag_load_balance.cc
  LIBRARIES oomph::solid oomph::constitutive oomph::navier_stokes oomph::meshes
            oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME mpi.multi_domain.turek_flag
  DEPENDS_ON turek_flag turek_flag_load_balance
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR} ${OOMPH_MPI_RUN_COMMAND}
  TEST_FILES validate.sh validata turek_flag_partition.dat load_balance.mcr
  LABELS mpi turek_flag)

# Don't run this test alongside any other tests
set_tests_properties(mpi.multi_domain.turek_flag
                     PROPERTIES PROCESSORS ${OOMPH_MPI_NUM_PROC})
if(APPLE)
  set_tests_properties(mpi.multi_domain.turek_flag PROPERTIES DISABLED FALSE)
endif()

message(VERBOSE "Leaving turek_flag subdirectory")
# ------------------------------------------------------------------------------
