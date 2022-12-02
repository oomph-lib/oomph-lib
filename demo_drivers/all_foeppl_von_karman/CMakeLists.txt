# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered all_foeppl_von_karman subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(all_foeppl_von_karman C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME axisym_displ_based_fvk
  SOURCES axisym_displ_based_fvk.cc
  LIBRARIES oomph::axisym_displ_based_foeppl_von_karman oomph::meshes
            oomph::generic)

oomph_add_executable(
  NAME axisym_fvk
  SOURCES axisym_fvk.cc
  LIBRARIES oomph::axisym_foeppl_von_karman oomph::meshes oomph::generic)

oomph_add_executable(
  NAME circular_disk
  SOURCES circular_disk.cc
  LIBRARIES oomph::foeppl_von_karman oomph::meshes oomph::generic)

oomph_add_executable(
  NAME displacement_based_circular_disk
  SOURCES displacement_based_circular_disk.cc
  LIBRARIES oomph::foeppl_von_karman oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME all_foeppl_von_karman.all_foeppl_von_karman
  TARGET_DEPENDENCIES axisym_displ_based_fvk axisym_fvk circular_disk
                      displacement_based_circular_disk
  EXTRA_REQUIRES
  LABELS all_foeppl_von_karman)

message(VERBOSE "Leaving all_foeppl_von_karman subdirectory")
# ------------------------------------------------------------------------------
