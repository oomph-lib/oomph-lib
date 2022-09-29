# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered two_d_adv_diff_SUPG subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(two_d_adv_diff_SUPG C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME two_d_adv_diff_SUPG
  SOURCES two_d_adv_diff_SUPG.cc
  LIBRARIES oomph::advection_diffusion
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME advection_diffusion.two_d_adv_diff_SUPG
  TARGET_DEPENDENCIES two_d_adv_diff_SUPG
  EXTRA_REQUIRES
  LABELS "advection_diffusion;two_d_adv_diff_SUPG")

message(VERBOSE "Leaving two_d_adv_diff_SUPG subdirectory")
# ------------------------------------------------------------------------------
