# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered two_d_unsteady_heat_t_adapt subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(two_d_unsteady_heat_t_adapt C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME two_d_unsteady_heat_t_adapt
  SOURCES two_d_unsteady_heat_t_adapt.cc
  LIBRARIES oomph::unsteady_heat)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME unsteady_heat.two_d_unsteady_heat_t_adapt
  TARGET_DEPENDENCIES two_d_unsteady_heat_t_adapt
  EXTRA_REQUIRES
  LABELS "unsteady_heat;two_d_unsteady_heat_t_adapt;adapt")

message(VERBOSE "Leaving two_d_unsteady_heat_t_adapt subdirectory")
# ------------------------------------------------------------------------------
