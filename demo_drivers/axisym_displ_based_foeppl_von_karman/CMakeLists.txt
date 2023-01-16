# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered axisym_displ_based_foeppl_von_karman subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(axisym_displ_based_foeppl_von_karman C CXX Fortran)
if(NOT OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME axisym_displ_based_fvk
  SOURCES axisym_displ_based_fvk.cc
  LIBRARIES oomph::axisym_displ_based_foeppl_von_karman oomph::meshes
            oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME
    axisym_displ_based_foeppl_von_karman.axisym_displ_based_foeppl_von_karman
  TARGET_DEPENDENCIES axisym_displ_based_fvk
  VALIDATE_SH_ARGS ${OOMPH_ROOT_DIR}
  EXTRA_REQUIRES
  LABELS axisym_displ_based_foeppl_von_karman)

message(VERBOSE "Leaving axisym_displ_based_foeppl_von_karman subdirectory")
# ------------------------------------------------------------------------------
