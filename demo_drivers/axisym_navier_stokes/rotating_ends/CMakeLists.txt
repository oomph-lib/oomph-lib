# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered rotating_ends subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(rotating_ends C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME rotating_ends
  SOURCES rotating_ends.cc
  LIBRARIES oomph::axisym_navier_stokes oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME axisym_navier_stokes.rotating_ends
  DEPENDS_ON rotating_ends
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving rotating_ends subdirectory")
# ------------------------------------------------------------------------------
