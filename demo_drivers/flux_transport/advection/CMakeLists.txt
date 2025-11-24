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
message(VERBOSE "Entered advection subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(advection C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()

include(CTest)

oomph_add_executable(
  NAME one_d_advection
  SOURCES one_d_advection.cc
  LIBRARIES oomph::flux_transport oomph::meshes oomph::generic)

oomph_add_executable(
  NAME two_d_advection
  SOURCES two_d_advection.cc
  LIBRARIES oomph::flux_transport oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME flux_transport.advection
  DEPENDS_ON one_d_advection two_d_advection
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving advection subdirectory")
# ------------------------------------------------------------------------------
