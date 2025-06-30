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
message(VERBOSE "Entered axisym_displ_based_foeppl_von_karman subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
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
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving axisym_displ_based_foeppl_von_karman subdirectory")
# ------------------------------------------------------------------------------
