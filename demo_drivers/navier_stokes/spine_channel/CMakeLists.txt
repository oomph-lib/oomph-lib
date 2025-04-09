# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered spine_channel subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(spine_channel C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME simple_spine_channel
  SOURCES simple_spine_channel.cc
  LIBRARIES oomph::navier_stokes oomph::meshes oomph::generic)

oomph_add_executable(
  NAME spine_channel
  SOURCES spine_channel.cc
  LIBRARIES oomph::navier_stokes oomph::meshes oomph::generic)

oomph_add_executable(
  NAME spine_channel2
  SOURCES spine_channel2.cc
  LIBRARIES oomph::navier_stokes oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME navier_stokes.spine_channel
  DEPENDS_ON simple_spine_channel spine_channel spine_channel2
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving spine_channel subdirectory")
# ------------------------------------------------------------------------------
