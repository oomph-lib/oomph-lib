# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered osc_ring subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(osc_ring C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME osc_ring_alg
  SOURCES osc_ring_alg.cc
  LIBRARIES oomph::navier_stokes oomph::meshes oomph::generic)

oomph_add_executable(
  NAME osc_ring_macro
  SOURCES osc_ring_macro.cc
  LIBRARIES oomph::navier_stokes oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME navier_stokes.osc_ring
  DEPENDS_ON osc_ring_alg osc_ring_macro
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving osc_ring subdirectory")
# ------------------------------------------------------------------------------
