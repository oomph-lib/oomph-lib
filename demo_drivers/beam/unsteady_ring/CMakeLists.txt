# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered unsteady_ring subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(unsteady_ring C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME lin_unsteady_ring
  SOURCES lin_unsteady_ring.cc
  LIBRARIES oomph::beam oomph::meshes oomph::generic)

oomph_add_executable(
  NAME unsteady_ring
  SOURCES unsteady_ring.cc
  LIBRARIES oomph::beam oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME beam.unsteady_ring
  DEPENDS_ON lin_unsteady_ring unsteady_ring
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving unsteady_ring subdirectory")
# ------------------------------------------------------------------------------
