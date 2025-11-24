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
message(VERBOSE "Entered fsi_collapsible_channel subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(fsi_collapsible_channel C CXX Fortran)

if(NOT oomphlib_FOUND)
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
  DEPENDS_ON fsi_collapsible_channel_algebraic_adapt_CR
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
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

# This test needs more than the default allotted time to complete
set_tests_properties(interaction.fsi_collapsible_channel PROPERTIES TIMEOUT
                                                                    3600)

message(VERBOSE "Leaving fsi_collapsible_channel subdirectory")
# ------------------------------------------------------------------------------
