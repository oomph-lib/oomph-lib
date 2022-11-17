# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered eluting_sphere subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(eluting_sphere C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME eluting_sphere
  SOURCES eluting_sphere.cc
  LIBRARIES oomph::spherical_advection_diffusion
  CXX_STANDARD 17)

oomph_add_executable(
  NAME eluting_sphere_adapt
  SOURCES eluting_sphere_adapt.cc
  LIBRARIES oomph::spherical_advection_diffusion
  CXX_STANDARD 17)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME spherical_advection_diffusion.eluting_sphere
  TARGET_DEPENDENCIES eluting_sphere eluting_sphere_adapt
  EXTRA_REQUIRES
  LABELS spherical_advection_diffusion eluting_sphere)

message(VERBOSE "Leaving eluting_sphere subdirectory")
# ------------------------------------------------------------------------------
