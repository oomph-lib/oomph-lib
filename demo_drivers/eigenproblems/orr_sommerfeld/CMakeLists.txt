# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered orr_sommerfeld subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(orr_sommerfeld C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME orr_sommerfeld
  SOURCES orr_sommerfeld.cc
  LIBRARIES oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME eigenproblems.orr_sommerfeld
  DEPENDS_ON orr_sommerfeld
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving orr_sommerfeld subdirectory")
# ------------------------------------------------------------------------------
