# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered curved_pipe subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(curved_pipe C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME curved_pipe
  SOURCES curved_pipe.cc
  LIBRARIES oomph::navier_stokes oomph::meshes oomph::generic)

oomph_add_executable(
  NAME helical_pipe
  SOURCES helical_pipe.cc
  LIBRARIES oomph::navier_stokes oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME navier_stokes.curved_pipe
  DEPENDS_ON curved_pipe helical_pipe
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving curved_pipe subdirectory")
# ------------------------------------------------------------------------------
