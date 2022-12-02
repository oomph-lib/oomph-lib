# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered pml_time_harmonic_linear_elasticity subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(pml_time_harmonic_linear_elasticity C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME time_harmonic_elasticity_driver
  SOURCES time_harmonic_elasticity_driver.cc
  LIBRARIES oomph::pml_time_harmonic_linear_elasticity oomph::meshes
            oomph::generic)

oomph_add_executable(
  NAME time_harmonic_elasticity_driver_source
  SOURCES time_harmonic_elasticity_driver_source.cc
  LIBRARIES oomph::pml_time_harmonic_linear_elasticity oomph::meshes
            oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME
    pml_time_harmonic_linear_elasticity.pml_time_harmonic_linear_elasticity
  TARGET_DEPENDENCIES time_harmonic_elasticity_driver
                      time_harmonic_elasticity_driver_source
  EXTRA_REQUIRES
  LABELS pml_time_harmonic_linear_elasticity)

message(VERBOSE "Leaving pml_time_harmonic_linear_elasticity subdirectory")
# ------------------------------------------------------------------------------
