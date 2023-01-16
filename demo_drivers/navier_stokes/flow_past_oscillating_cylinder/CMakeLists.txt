# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered flow_past_oscillating_cylinder subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(flow_past_oscillating_cylinder C CXX Fortran)
if(NOT OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME flow_past_oscillating_cylinder
  SOURCES flow_past_oscillating_cylinder.cc
  LIBRARIES oomph::navier_stokes oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME navier_stokes.flow_past_oscillating_cylinder
  TARGET_DEPENDENCIES flow_past_oscillating_cylinder
  VALIDATE_SH_ARGS ${OOMPH_ROOT_DIR}
  EXTRA_REQUIRES
  LABELS navier_stokes flow_past_oscillating_cylinder)

message(VERBOSE "Leaving flow_past_oscillating_cylinder subdirectory")
# ------------------------------------------------------------------------------
