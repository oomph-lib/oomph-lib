# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered two_d_act_inhibit subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(two_d_act_inhibit C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME two_d_act_inhibit
  SOURCES two_d_act_inhibit.cc
  LIBRARIES oomph::advection_diffusion_reaction
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME reaction_diffusion.two_d_act_inhibit
  TARGET_DEPENDENCIES two_d_act_inhibit
  EXTRA_REQUIRES
  LABELS reaction_diffusion two_d_act_inhibit)

message(VERBOSE "Leaving two_d_act_inhibit subdirectory")
# ------------------------------------------------------------------------------
