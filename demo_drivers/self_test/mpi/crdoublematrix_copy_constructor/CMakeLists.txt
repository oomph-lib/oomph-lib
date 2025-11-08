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
message(VERBOSE "Entered crdoublematrix_copy_constructor subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(crdoublematrix_copy_constructor C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../../install")
endif()

include(CTest)

oomph_add_executable(
  NAME crdoublematrix_copy_constructor
  SOURCES crdoublematrix_copy_constructor.cc
  LIBRARIES oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

include(ProcessorCount)
ProcessorCount(NUM_PROCESSORS)
if(NUM_PROCESSORS GREATER_EQUAL 4)
  oomph_add_test(
    TEST_NAME self_test.mpi.crdoublematrix_copy_constructor
    DEPENDS_ON crdoublematrix_copy_constructor
    COMMAND ./validate.sh ${OOMPH_ROOT_DIR} ${OOMPH_MPI_VARIABLENP_RUN_COMMAND}
    TEST_FILES validate.sh validata)

  # Don't run this test alongside any other tests
  set_tests_properties(self_test.mpi.crdoublematrix_copy_constructor
                       PROPERTIES PROCESSORS 4)
else()
  message(
    STATUS
      "Skipping test 'self_test.mpi.crdoublematrix_copy_constructor' as you don't have enough processors available."
  )
endif()

message(VERBOSE "Leaving crdoublematrix_copy_constructor subdirectory")
# ------------------------------------------------------------------------------
