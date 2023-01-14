# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered pseudo_solid_collapsible_tube subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(pseudo_solid_collapsible_tube C CXX Fortran)
if(NOT OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME pseudo_solid_collapsible_tube
  SOURCES pseudo_solid_collapsible_tube.cc
  LIBRARIES oomph::multi_physics oomph::navier_stokes oomph::solid
            oomph::constitutive oomph::meshes oomph::generic)

oomph_add_executable(
  NAME pseudo_solid_collapsible_tube_work
  SOURCES pseudo_solid_collapsible_tube_work.cc
  LIBRARIES oomph::multi_physics oomph::navier_stokes oomph::solid
            oomph::constitutive oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME interaction.pseudo_solid_collapsible_tube
  TARGET_DEPENDENCIES pseudo_solid_collapsible_tube
                      pseudo_solid_collapsible_tube_work
  VALIDATE_SH_ARGS ${OOMPH_ROOT_DIR}
  EXTRA_REQUIRES
  LABELS interaction pseudo_solid_collapsible_tube)

message(VERBOSE "Leaving pseudo_solid_collapsible_tube subdirectory")
# ------------------------------------------------------------------------------
