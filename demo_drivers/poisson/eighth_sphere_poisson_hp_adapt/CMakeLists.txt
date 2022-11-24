# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered eighth_sphere_poisson_hp_adapt subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(eighth_sphere_poisson_hp_adapt C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME eighth_sphere_poisson_hp_adapt
  SOURCES eighth_sphere_poisson_hp_adapt.cc
  LIBRARIES oomph::poisson)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME poisson.eighth_sphere_poisson_hp_adapt
  TARGET_DEPENDENCIES eighth_sphere_poisson_hp_adapt
  EXTRA_REQUIRES
  LABELS "poisson;eighth_sphere_poisson_hp_adapt;adapt")

message(VERBOSE "Leaving eighth_sphere_poisson_hp_adapt subdirectory")
# ------------------------------------------------------------------------------
