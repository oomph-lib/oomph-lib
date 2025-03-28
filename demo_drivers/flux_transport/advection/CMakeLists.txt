# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered advection subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(advection C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME one_d_advection
  SOURCES one_d_advection.cc
  LIBRARIES oomph::flux_transport oomph::meshes oomph::generic)

oomph_add_executable(
  NAME two_d_advection
  SOURCES two_d_advection.cc
  LIBRARIES oomph::flux_transport oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME flux_transport.advection
  DEPENDS_ON one_d_advection two_d_advection
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving advection subdirectory")
# ------------------------------------------------------------------------------
