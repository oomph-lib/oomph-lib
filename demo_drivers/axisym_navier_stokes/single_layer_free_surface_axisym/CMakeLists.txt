# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered single_layer_free_surface_axisym subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(single_layer_free_surface_axisym C CXX Fortran)
if(NOT OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME single_layer_free_surface_axisym
  SOURCES single_layer_free_surface_axisym.cc
  LIBRARIES oomph::axisym_navier_stokes oomph::navier_stokes
            oomph::fluid_interface oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME axisym_navier_stokes.single_layer_free_surface_axisym
  TARGET_DEPENDENCIES single_layer_free_surface_axisym
  EXTRA_REQUIRES
  LABELS axisym_navier_stokes single_layer_free_surface_axisym)

message(VERBOSE "Leaving single_layer_free_surface_axisym subdirectory")
# ------------------------------------------------------------------------------
