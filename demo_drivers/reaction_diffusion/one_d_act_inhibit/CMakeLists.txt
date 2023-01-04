# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered one_d_act_inhibit subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(one_d_act_inhibit C CXX Fortran)
if(NOT OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME one_d_act_inhibit
  SOURCES one_d_act_inhibit.cc
  LIBRARIES oomph::advection_diffusion_reaction oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME reaction_diffusion.one_d_act_inhibit
  TARGET_DEPENDENCIES one_d_act_inhibit
  EXTRA_REQUIRES
  LABELS reaction_diffusion one_d_act_inhibit)

message(VERBOSE "Leaving one_d_act_inhibit subdirectory")
# ------------------------------------------------------------------------------
