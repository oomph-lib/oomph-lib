# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered euler subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(euler C CXX Fortran)
if(NOT OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME couette
  SOURCES couette.cc
  LIBRARIES oomph::flux_transport oomph::meshes oomph::generic)

oomph_add_executable(
  NAME ff_step
  SOURCES ff_step.cc
  LIBRARIES oomph::flux_transport oomph::meshes oomph::generic)

oomph_add_executable(
  NAME one_d_euler
  SOURCES one_d_euler.cc
  LIBRARIES oomph::flux_transport oomph::meshes oomph::generic)

oomph_add_executable(
  NAME ramp
  SOURCES ramp.cc
  LIBRARIES oomph::flux_transport oomph::meshes oomph::generic)

oomph_add_executable(
  NAME two_d_euler
  SOURCES two_d_euler.cc
  LIBRARIES oomph::flux_transport oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME flux_transport.euler
  TARGET_DEPENDENCIES couette ff_step one_d_euler ramp two_d_euler
  VALIDATE_SH_ARGS ${OOMPH_ROOT_DIR}
  EXTRA_REQUIRES
  LABELS flux_transport euler)

message(VERBOSE "Leaving euler subdirectory")
# ------------------------------------------------------------------------------
