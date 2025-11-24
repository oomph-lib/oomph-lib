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
message(VERBOSE "Entered two_layer_soluble_surfactant subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(two_layer_soluble_surfactant C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()

include(CTest)

oomph_add_executable(
  NAME refineable_two_layer_soluble_surfactant
  SOURCES double_buoyant_navier_stokes_elements.h
          refineable_two_layer_soluble_surfactant.cc
  LIBRARIES oomph::advection_diffusion_reaction
            oomph::navier_stokes
            oomph::fluid_interface
            oomph::solid
            oomph::constitutive
            oomph::meshes
            oomph::generic)

oomph_add_executable(
  NAME two_layer_soluble_surfactant
  SOURCES soluble_surfactant_transport_equations.h
          double_buoyant_navier_stokes_elements.h
          two_layer_soluble_surfactant.cc
          soluble_surfactant_transport_equations.h
  LIBRARIES oomph::advection_diffusion_reaction oomph::navier_stokes
            oomph::fluid_interface oomph::meshes oomph::generic)

oomph_add_executable(
  NAME refineable_unstructured_two_layer_ss
  SOURCES double_buoyant_navier_stokes_elements.h
          refineable_unstructured_two_layer_ss.cc
  LIBRARIES oomph::advection_diffusion_reaction
            oomph::navier_stokes
            oomph::fluid_interface
            oomph::solid
            oomph::constitutive
            oomph::meshes
            oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME multi_physics.two_layer_soluble_surfactant
  DEPENDS_ON refineable_two_layer_soluble_surfactant
             two_layer_soluble_surfactant refineable_unstructured_two_layer_ss
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving two_layer_soluble_surfactant subdirectory")
# ------------------------------------------------------------------------------
