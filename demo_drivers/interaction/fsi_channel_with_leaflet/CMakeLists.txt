# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered fsi_channel_with_leaflet subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(fsi_channel_with_leaflet C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME fsi_channel_with_leaflet
  SOURCES fsi_channel_with_leaflet.cc
  LIBRARIES oomph::multi_physics oomph::beam oomph::navier_stokes oomph::meshes
            oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME interaction.fsi_channel_with_leaflet
  DEPENDS_ON fsi_channel_with_leaflet
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving fsi_channel_with_leaflet subdirectory")
# ------------------------------------------------------------------------------
