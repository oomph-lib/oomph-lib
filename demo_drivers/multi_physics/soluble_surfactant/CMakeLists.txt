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
message(VERBOSE "Entered soluble_surfactant subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(soluble_surfactant C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()

include(CTest)

oomph_add_executable(
  NAME soluble_surfactant
  SOURCES double_buoyant_navier_stokes_elements.h soluble_surfactant.cc
  LIBRARIES oomph::advection_diffusion_reaction oomph::navier_stokes
            oomph::fluid_interface oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME multi_physics.soluble_surfactant
  DEPENDS_ON soluble_surfactant
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving soluble_surfactant subdirectory")
# ------------------------------------------------------------------------------
