# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered two_d_adv_diff_adapt subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(two_d_adv_diff_adapt C CXX Fortran)
if(NOT OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME two_d_adv_diff_adapt
  SOURCES two_d_adv_diff_adapt.cc
  LIBRARIES oomph::advection_diffusion oomph::meshes oomph::generic)

oomph_add_executable(
  NAME two_d_adv_diff_adapt2
  SOURCES two_d_adv_diff_adapt2.cc
  LIBRARIES oomph::advection_diffusion oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME advection_diffusion.two_d_adv_diff_adapt
  TARGET_DEPENDENCIES two_d_adv_diff_adapt two_d_adv_diff_adapt2
  EXTRA_REQUIRES
  LABELS advection_diffusion two_d_adv_diff_adapt adapt)

message(VERBOSE "Leaving two_d_adv_diff_adapt subdirectory")
# ------------------------------------------------------------------------------
