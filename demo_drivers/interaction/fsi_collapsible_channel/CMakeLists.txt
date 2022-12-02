# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered fsi_collapsible_channel subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(fsi_collapsible_channel C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME fsi_collapsible_channel_algebraic_adapt_CR
  SOURCES fsi_collapsible_channel_adapt.cc
  LIBRARIES oomph::beam oomph::navier_stokes oomph::meshes oomph::generic)

oomph_add_executable(
  NAME fsi_collapsible_channel_algebraic_adapt_TH
  SOURCES fsi_collapsible_channel_adapt.cc
  LIBRARIES oomph::beam oomph::navier_stokes oomph::meshes oomph::generic
  CXX_DEFINITIONS TAYLOR_HOOD)

oomph_add_executable(
  NAME fsi_collapsible_channel_algebraic_CR
  SOURCES fsi_collapsible_channel.cc
  LIBRARIES oomph::beam oomph::navier_stokes oomph::meshes oomph::generic)

oomph_add_executable(
  NAME fsi_collapsible_channel_algebraic_TH
  SOURCES fsi_collapsible_channel.cc
  LIBRARIES oomph::beam oomph::navier_stokes oomph::meshes oomph::generic
  CXX_DEFINITIONS TAYLOR_HOOD)

oomph_add_executable(
  NAME fsi_collapsible_channel_macro_adapt_CR
  SOURCES fsi_collapsible_channel_adapt.cc
  LIBRARIES oomph::beam oomph::navier_stokes oomph::meshes oomph::generic
  CXX_DEFINITIONS MACRO_ELEMENT_NODE_UPDATE)

oomph_add_executable(
  NAME fsi_collapsible_channel_macro_adapt_TH
  SOURCES fsi_collapsible_channel_adapt.cc
  LIBRARIES oomph::beam oomph::navier_stokes oomph::meshes oomph::generic
  CXX_DEFINITIONS MACRO_ELEMENT_NODE_UPDATE TAYLOR_HOOD)

oomph_add_executable(
  NAME fsi_collapsible_channel_macro_CR
  SOURCES fsi_collapsible_channel.cc
  LIBRARIES oomph::beam oomph::navier_stokes oomph::meshes oomph::generic
  CXX_DEFINITIONS MACRO_ELEMENT_NODE_UPDATE)

oomph_add_executable(
  NAME fsi_collapsible_channel_macro_TH
  SOURCES fsi_collapsible_channel.cc
  LIBRARIES oomph::beam oomph::navier_stokes oomph::meshes oomph::generic
  CXX_DEFINITIONS MACRO_ELEMENT_NODE_UPDATE TAYLOR_HOOD)

oomph_add_executable(
  NAME fsi_pseudo_solid_collapsible_channel_adapt_CR
  SOURCES fsi_pseudo_solid_collapsible_channel_adapt.cc
  LIBRARIES oomph::beam oomph::constitutive oomph::solid oomph::navier_stokes
            oomph::meshes oomph::generic)

oomph_add_executable(
  NAME fsi_pseudo_solid_collapsible_channel_adapt_TH
  SOURCES fsi_pseudo_solid_collapsible_channel_adapt.cc
  LIBRARIES oomph::beam oomph::constitutive oomph::solid oomph::navier_stokes
            oomph::meshes oomph::generic
  CXX_DEFINITIONS TAYLOR_HOOD)

oomph_add_executable(
  NAME fsi_pseudo_solid_collapsible_channel_CR
  SOURCES fsi_pseudo_solid_collapsible_channel.cc
  LIBRARIES oomph::beam oomph::constitutive oomph::solid oomph::navier_stokes
            oomph::meshes oomph::generic)

oomph_add_executable(
  NAME fsi_pseudo_solid_collapsible_channel_TH
  SOURCES fsi_pseudo_solid_collapsible_channel.cc
  LIBRARIES oomph::beam oomph::constitutive oomph::solid oomph::navier_stokes
            oomph::meshes oomph::generic
  CXX_DEFINITIONS TAYLOR_HOOD)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME interaction.fsi_collapsible_channel
  TARGET_DEPENDENCIES
    fsi_collapsible_channel_algebraic_adapt_CR
    fsi_collapsible_channel_algebraic_adapt_TH
    fsi_collapsible_channel_algebraic_CR
    fsi_collapsible_channel_algebraic_TH
    fsi_collapsible_channel_macro_adapt_CR
    fsi_collapsible_channel_macro_adapt_TH
    fsi_collapsible_channel_macro_CR
    fsi_collapsible_channel_macro_TH
    fsi_pseudo_solid_collapsible_channel_adapt_CR
    fsi_pseudo_solid_collapsible_channel_adapt_TH
    fsi_pseudo_solid_collapsible_channel_CR
    fsi_pseudo_solid_collapsible_channel_TH
  EXTRA_REQUIRES
  LABELS interaction fsi_collapsible_channel)

# This test needs more than the default allotted time to complete
set_tests_properties(interaction.fsi_collapsible_channel PROPERTIES TIMEOUT
                                                                    10000)

message(VERBOSE "Leaving fsi_collapsible_channel subdirectory")
# ------------------------------------------------------------------------------
