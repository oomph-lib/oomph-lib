# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered hp_adaptive_driven_cavity subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(hp_adaptive_driven_cavity C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME circular_driven_cavity_hp_adapt
  SOURCES circular_driven_cavity_hp_adapt.cc
  LIBRARIES oomph::navier_stokes oomph::meshes oomph::generic)

oomph_add_executable(
  NAME hp_adaptive_driven_cavity
  SOURCES hp_adaptive_driven_cavity.cc
  LIBRARIES oomph::navier_stokes oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME navier_stokes.hp_adaptive_driven_cavity
  DEPENDS_ON circular_driven_cavity_hp_adapt hp_adaptive_driven_cavity
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata
  LABELS navier_stokes hp_adaptive_driven_cavity adapt)

message(VERBOSE "Leaving hp_adaptive_driven_cavity subdirectory")
# ------------------------------------------------------------------------------
