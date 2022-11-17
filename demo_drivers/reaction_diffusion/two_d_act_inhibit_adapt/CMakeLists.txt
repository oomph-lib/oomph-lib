# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered two_d_act_inhibit_adapt subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(two_d_act_inhibit_adapt C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME two_d_act_inhibit_adapt
  SOURCES two_d_act_inhibit_adapt.cc
  LIBRARIES oomph::advection_diffusion_reaction
  CXX_STANDARD 17)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME reaction_diffusion.two_d_act_inhibit_adapt
  TARGET_DEPENDENCIES two_d_act_inhibit_adapt
  EXTRA_REQUIRES
  LABELS "reaction_diffusion;two_d_act_inhibit_adapt;adapt")

message(VERBOSE "Leaving two_d_act_inhibit_adapt subdirectory")
# ------------------------------------------------------------------------------
