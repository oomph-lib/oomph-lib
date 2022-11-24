# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered marangoni_convection subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(marangoni_convection C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME marangoni_convection
  SOURCES marangoni_convection.cc
  LIBRARIES oomph::advection_diffusion oomph::navier_stokes
            oomph::fluid_interface oomph::multi_physics)

oomph_add_executable(
  NAME marangoni_convection_box
  SOURCES marangoni_convection_box.cc
  LIBRARIES oomph::advection_diffusion oomph::navier_stokes
            oomph::multi_physics oomph::fluid_interface)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME multi_physics.marangoni_convection
  TARGET_DEPENDENCIES marangoni_convection marangoni_convection_box
  EXTRA_REQUIRES
  LABELS multi_physics marangoni_convection)

message(VERBOSE "Leaving marangoni_convection subdirectory")
# ------------------------------------------------------------------------------
