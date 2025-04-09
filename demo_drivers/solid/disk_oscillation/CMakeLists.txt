# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered disk_oscillation subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(disk_oscillation C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME disk_oscillation
  SOURCES disk_oscillation.cc
  LIBRARIES oomph::solid oomph::constitutive oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME solid.disk_oscillation
  DEPENDS_ON disk_oscillation
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving disk_oscillation subdirectory")
# ------------------------------------------------------------------------------
