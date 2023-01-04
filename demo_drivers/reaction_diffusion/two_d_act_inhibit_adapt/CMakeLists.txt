# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered two_d_act_inhibit_adapt subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(two_d_act_inhibit_adapt C CXX Fortran)
if(NOT OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME two_d_act_inhibit_adapt
  SOURCES two_d_act_inhibit_adapt.cc
  LIBRARIES oomph::advection_diffusion_reaction oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME reaction_diffusion.two_d_act_inhibit_adapt
  TARGET_DEPENDENCIES two_d_act_inhibit_adapt
  EXTRA_REQUIRES
  LABELS reaction_diffusion two_d_act_inhibit_adapt adapt)

message(VERBOSE "Leaving two_d_act_inhibit_adapt subdirectory")
# ------------------------------------------------------------------------------
