# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered osc_ellipse subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(osc_ellipse C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME osc_quarter_ellipse
  SOURCES osc_quarter_ellipse.cc
  LIBRARIES oomph::navier_stokes oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME navier_stokes.osc_ellipse
  DEPENDS_ON osc_quarter_ellipse
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving osc_ellipse subdirectory")
# ------------------------------------------------------------------------------
