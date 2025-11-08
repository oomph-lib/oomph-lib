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
message(VERBOSE "Entered pseudo_solid_fsi_channel_with_leaflet subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(pseudo_solid_fsi_channel_with_leaflet C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()

include(CTest)

oomph_add_executable(
  NAME fsi_channel_with_leaflet_precond
  SOURCES fsi_channel_with_leaflet_precond.cc
  LIBRARIES oomph::multi_physics
            oomph::beam
            oomph::navier_stokes
            oomph::constitutive
            oomph::solid
            oomph::meshes
            oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME interaction.pseudo_solid_fsi_channel_with_leaflet
  DEPENDS_ON fsi_channel_with_leaflet_precond
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving pseudo_solid_fsi_channel_with_leaflet subdirectory")
# ------------------------------------------------------------------------------
