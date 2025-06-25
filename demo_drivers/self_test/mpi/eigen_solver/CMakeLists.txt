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
message(VERBOSE "Entered eigen_solver subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(eigen_solver C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../../install")
endif()

include(CTest)

if(OOMPH_HAS_TRILINOS)
  oomph_add_executable(
    NAME eigen_solver_test
    SOURCES eigen_solver_test.cc
    LIBRARIES oomph::meshes oomph::generic)
endif()

# ---------------------------------- TESTING ----------------------------------

if(OOMPH_HAS_TRILINOS)
  oomph_add_test(
    TEST_NAME self_test.mpi.eigen_solver
    DEPENDS_ON eigen_solver
    COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
    TEST_FILES validate.sh validata random_test_matrix.dat)
endif()

message(VERBOSE "Leaving eigen_solver subdirectory")
# ------------------------------------------------------------------------------
