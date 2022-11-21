# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered collapsible_channel subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(collapsible_channel C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME collapsible_channel
  SOURCES collapsible_channel.cc
  LIBRARIES oomph::navier_stokes
  CXX_STANDARD 17)

oomph_add_executable(
  NAME collapsible_channel_adaptive_algebraic
  SOURCES collapsible_channel_algebraic.cc my_alg_channel_mesh.h
  LIBRARIES oomph::navier_stokes
  CXX_DEFINITIONS ADAPTIVE
  CXX_STANDARD 17)

oomph_add_executable(
  NAME collapsible_channel_adaptive_algebraic_bl_squash
  SOURCES collapsible_channel_algebraic.cc my_alg_channel_mesh.h
  LIBRARIES oomph::navier_stokes
  CXX_DEFINITIONS USE_BL_SQUASH_FCT ADAPTIVE
  CXX_STANDARD 17)

oomph_add_executable(
  NAME collapsible_channel_algebraic
  SOURCES collapsible_channel_algebraic.cc my_alg_channel_mesh.h
  LIBRARIES oomph::navier_stokes
  CXX_STANDARD 17)

oomph_add_executable(
  NAME collapsible_channel_algebraic_bl_squash
  SOURCES collapsible_channel_algebraic.cc my_alg_channel_mesh.h
  LIBRARIES oomph::navier_stokes
  CXX_DEFINITIONS USE_BL_SQUASH_FCT
  CXX_STANDARD 17)

oomph_add_executable(
  NAME collapsible_channel_bl_squash
  SOURCES collapsible_channel.cc
  LIBRARIES oomph::navier_stokes
  CXX_DEFINITIONS USE_BL_SQUASH_FCT
  CXX_STANDARD 17)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME navier_stokes.collapsible_channel
  TARGET_DEPENDENCIES
    collapsible_channel
    collapsible_channel_adaptive_algebraic
    collapsible_channel_adaptive_algebraic_bl_squash
    collapsible_channel_algebraic
    collapsible_channel_algebraic_bl_squash
    collapsible_channel_bl_squash
  EXTRA_REQUIRES
  LABELS navier_stokes collapsible_channel)

message(VERBOSE "Leaving collapsible_channel subdirectory")
# ------------------------------------------------------------------------------
