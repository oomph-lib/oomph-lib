# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered rayleigh_traction_channel subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(rayleigh_traction_channel C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME rayleigh_traction_channel
  SOURCES rayleigh_traction_channel.cc
  LIBRARIES oomph::navier_stokes oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME navier_stokes.rayleigh_traction_channel
  TARGET_DEPENDENCIES rayleigh_traction_channel
  VALIDATE_SH_ARGS ${OOMPH_ROOT_DIR}
  EXTRA_REQUIRES
  LABELS navier_stokes rayleigh_traction_channel)

message(VERBOSE "Leaving rayleigh_traction_channel subdirectory")
# ------------------------------------------------------------------------------
