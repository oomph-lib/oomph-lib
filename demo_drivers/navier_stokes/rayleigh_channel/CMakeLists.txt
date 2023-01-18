# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered rayleigh_channel subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(rayleigh_channel C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME rayleigh_channel
  SOURCES rayleigh_channel.cc
  LIBRARIES oomph::navier_stokes oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME navier_stokes.rayleigh_channel
  DEPENDS_ON rayleigh_channel
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata
  LABELS navier_stokes rayleigh_channel)

message(VERBOSE "Leaving rayleigh_channel subdirectory")
# ------------------------------------------------------------------------------
