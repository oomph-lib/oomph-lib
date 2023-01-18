# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered axisym_displ_based_foeppl_von_karman subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(axisym_displ_based_foeppl_von_karman C CXX Fortran)
if(NOT oomphlib_FOUND)
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
  DEPENDS_ON axisym_displ_based_fvk
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata
  LABELS axisym_displ_based_foeppl_von_karman)

message(VERBOSE "Leaving axisym_displ_based_foeppl_von_karman subdirectory")
# ------------------------------------------------------------------------------
