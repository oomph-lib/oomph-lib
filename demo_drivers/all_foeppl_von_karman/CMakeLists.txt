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
message(VERBOSE "Entered all_foeppl_von_karman subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(all_foeppl_von_karman C CXX Fortran)

if(NOT oomphlib_FOUND)
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
  DEPENDS_ON axisym_displ_based_fvk axisym_fvk circular_disk
             displacement_based_circular_disk
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving all_foeppl_von_karman subdirectory")
# ------------------------------------------------------------------------------
