# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered navier_stokes_with_singularity subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(navier_stokes_with_singularity C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME adaptive_driven_cavity
  SOURCES adaptive_driven_cavity.cc
  LIBRARIES oomph::navier_stokes oomph::meshes oomph::generic)

oomph_add_executable(
  NAME circular_couette
  SOURCES driven_cavity.cc navier_stokes_elements_with_singularity.h
  LIBRARIES oomph::navier_stokes oomph::meshes oomph::generic
  CXX_DEFINITIONS COUETTE)

oomph_add_executable(
  NAME circular_couette_without_pseudo_singularity_subtracted_off
  SOURCES driven_cavity.cc navier_stokes_elements_with_singularity.h
          finite_re_perturbation.h
  LIBRARIES oomph::navier_stokes oomph::meshes oomph::generic
  CXX_DEFINITIONS DONT_USE_SINGULARITY COUETTE)

oomph_add_executable(
  NAME driven_cavity
  SOURCES driven_cavity.cc navier_stokes_elements_with_singularity.h
          finite_re_perturbation.h
  LIBRARIES oomph::navier_stokes oomph::meshes oomph::generic)

oomph_add_executable(
  NAME driven_cavity_without_singularity_subtracted_off
  SOURCES driven_cavity.cc navier_stokes_elements_with_singularity.h
          finite_re_perturbation.h
  LIBRARIES oomph::navier_stokes oomph::meshes oomph::generic
  CXX_DEFINITIONS DONT_USE_SINGULARITY)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME navier_stokes.navier_stokes_with_singularity
  TARGET_DEPENDENCIES
    adaptive_driven_cavity circular_couette
    circular_couette_without_pseudo_singularity_subtracted_off driven_cavity
    driven_cavity_without_singularity_subtracted_off
  VALIDATE_SH_ARGS ${OOMPH_ROOT_DIR}
  EXTRA_REQUIRES
  LABELS navier_stokes navier_stokes_with_singularity)

message(VERBOSE "Leaving navier_stokes_with_singularity subdirectory")
# ------------------------------------------------------------------------------
