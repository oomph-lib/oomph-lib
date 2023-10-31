# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered steady_ring subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(steady_ring C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME steady_ring
  SOURCES steady_ring.cc
  LIBRARIES oomph::beam oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME beam.steady_ring
  DEPENDS_ON steady_ring
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata
  LABELS beam steady_ring)

message(VERBOSE "Leaving steady_ring subdirectory")
# ------------------------------------------------------------------------------
