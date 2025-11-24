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
message(VERBOSE "Entered vector_split subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(vector_split C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../../install")
endif()

include(CTest)

oomph_add_executable(
  NAME vector_split
  SOURCES vector_split.cc
  LIBRARIES oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

include(ProcessorCount)
ProcessorCount(NUM_PROCESSORS)
if(NUM_PROCESSORS GREATER_EQUAL 4)
  oomph_add_test(
    TEST_NAME self_test.mpi.vector_split
    DEPENDS_ON vector_split
    COMMAND ./validate.sh ${OOMPH_ROOT_DIR} ${OOMPH_MPI_VARIABLENP_RUN_COMMAND}
    TEST_FILES validate.sh validata)

  # Don't run this test alongside any other tests
  set_tests_properties(self_test.mpi.vector_split PROPERTIES PROCESSORS 4)
else()
  message(
    STATUS
      "Skipping test 'self_test.mpi.vector_split' as you don't have enough processors available."
  )
endif()

message(VERBOSE "Leaving vector_split subdirectory")
# ------------------------------------------------------------------------------
