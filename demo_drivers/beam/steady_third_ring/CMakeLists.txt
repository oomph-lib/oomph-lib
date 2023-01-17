# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered steady_third_ring subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(steady_third_ring C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME steady_third_ring
  SOURCES steady_third_ring.cc
  LIBRARIES oomph::beam oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME beam.steady_third_ring
  TARGET_DEPENDENCIES steady_third_ring
  VALIDATE_SH_ARGS ${OOMPH_ROOT_DIR}
  EXTRA_REQUIRES
  LABELS beam steady_third_ring)

message(VERBOSE "Leaving steady_third_ring subdirectory")
# ------------------------------------------------------------------------------
