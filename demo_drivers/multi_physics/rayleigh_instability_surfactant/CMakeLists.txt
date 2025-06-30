# ------------------------------------------------------------------------------
# OOMPH auto-inserted: warn once if user configures inside an existing build.
if(EXISTS "${CMAKE_CURRENT_BINARY_DIR}/CMakeCache.txt"
   AND CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  message(
    WARNING
      "You're re-configuring an existing build tree. For a clean configuration, "
      "delete the whole build folder and start with a fresh directory.")
endif()
# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered rayleigh_instability_surfactant subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(rayleigh_instability_surfactant C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()

include(CTest)

oomph_add_executable(
  NAME 3d_rayleigh_instability_surfactant
  SOURCES 3d_rayleigh_instability_surfactant.cc
  LIBRARIES oomph::navier_stokes oomph::fluid_interface oomph::meshes
            oomph::generic)

oomph_add_executable(
  NAME rayleigh_instability_insoluble_surfactant
  SOURCES rayleigh_instability_insoluble_surfactant.cc
          axisymmetric_advection_navier_stokes_elements.h
  LIBRARIES oomph::axisym_navier_stokes oomph::fluid_interface oomph::meshes
            oomph::generic)

oomph_add_executable(
  NAME rayleigh_instability_soluble_surfactant
  SOURCES rayleigh_instability_soluble_surfactant.cc
  LIBRARIES oomph::axisym_navier_stokes
            oomph::fluid_interface
            oomph::solid
            oomph::constitutive
            oomph::axisym_advection_diffusion
            oomph::meshes
            oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME multi_physics.rayleigh_instability_surfactant
  DEPENDS_ON 3d_rayleigh_instability_surfactant
             rayleigh_instability_insoluble_surfactant
             rayleigh_instability_soluble_surfactant
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving rayleigh_instability_surfactant subdirectory")
# ------------------------------------------------------------------------------
