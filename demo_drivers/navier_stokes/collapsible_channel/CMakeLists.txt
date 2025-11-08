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
message(VERBOSE "Entered collapsible_channel subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(collapsible_channel C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()

include(CTest)

oomph_add_executable(
  NAME collapsible_channel
  SOURCES collapsible_channel.cc
  LIBRARIES oomph::navier_stokes oomph::meshes oomph::generic)

oomph_add_executable(
  NAME collapsible_channel_adaptive_algebraic
  SOURCES collapsible_channel_algebraic.cc my_alg_channel_mesh.h
  LIBRARIES oomph::navier_stokes oomph::meshes oomph::generic
  CXX_DEFINITIONS ADAPTIVE)

oomph_add_executable(
  NAME collapsible_channel_adaptive_algebraic_bl_squash
  SOURCES collapsible_channel_algebraic.cc my_alg_channel_mesh.h
  LIBRARIES oomph::navier_stokes oomph::meshes oomph::generic
  CXX_DEFINITIONS USE_BL_SQUASH_FCT ADAPTIVE)

oomph_add_executable(
  NAME collapsible_channel_algebraic
  SOURCES collapsible_channel_algebraic.cc my_alg_channel_mesh.h
  LIBRARIES oomph::navier_stokes oomph::meshes oomph::generic)

oomph_add_executable(
  NAME collapsible_channel_algebraic_bl_squash
  SOURCES collapsible_channel_algebraic.cc my_alg_channel_mesh.h
  LIBRARIES oomph::navier_stokes oomph::meshes oomph::generic
  CXX_DEFINITIONS USE_BL_SQUASH_FCT)

oomph_add_executable(
  NAME collapsible_channel_bl_squash
  SOURCES collapsible_channel.cc
  LIBRARIES oomph::navier_stokes oomph::meshes oomph::generic
  CXX_DEFINITIONS USE_BL_SQUASH_FCT)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME navier_stokes.collapsible_channel
  DEPENDS_ON collapsible_channel
             collapsible_channel_adaptive_algebraic
             collapsible_channel_adaptive_algebraic_bl_squash
             collapsible_channel_algebraic
             collapsible_channel_algebraic_bl_squash
             collapsible_channel_bl_squash
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving collapsible_channel subdirectory")
# ------------------------------------------------------------------------------
