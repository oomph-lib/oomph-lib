# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered steady_axisym_advection_diffusion subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(steady_axisym_advection_diffusion C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME steady_axisym_advection_diffusion
  SOURCES steady_axisym_advection_diffusion.cc
  LIBRARIES oomph::steady_axisym_advection_diffusion
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME steady_axisym_advection_diffusion.steady_axisym_advection_diffusion
  TARGET_DEPENDENCIES steady_axisym_advection_diffusion
  EXTRA_REQUIRES
  LABELS "steady_axisym_advection_diffusion;steady_axisym_advection_diffusion")

message(VERBOSE "Leaving steady_axisym_advection_diffusion subdirectory")
# ------------------------------------------------------------------------------
