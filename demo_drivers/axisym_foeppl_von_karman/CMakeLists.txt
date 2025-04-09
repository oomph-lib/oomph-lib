# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered axisym_foeppl_von_karman subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(axisym_foeppl_von_karman C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME axisym_fvk
  SOURCES axisym_fvk.cc
  LIBRARIES oomph::axisym_foeppl_von_karman oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME axisym_foeppl_von_karman.axisym_foeppl_von_karman
  DEPENDS_ON axisym_fvk
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving axisym_foeppl_von_karman subdirectory")
# ------------------------------------------------------------------------------
