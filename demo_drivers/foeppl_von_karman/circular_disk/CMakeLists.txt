# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered circular_disk subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(circular_disk C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME circular_disk
  SOURCES circular_disk.cc
  LIBRARIES oomph::foeppl_von_karman oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME foeppl_von_karman.circular_disk
  DEPENDS_ON circular_disk
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving circular_disk subdirectory")
# ------------------------------------------------------------------------------
