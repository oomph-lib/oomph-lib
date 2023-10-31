# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered spherical_couette subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(spherical_couette C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME spherical_couette
  SOURCES spherical_couette.cc
  LIBRARIES oomph::spherical_navier_stokes oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME spherical_navier_stokes.spherical_couette
  DEPENDS_ON spherical_couette
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata
  LABELS spherical_navier_stokes spherical_couette)

message(VERBOSE "Leaving spherical_couette subdirectory")
# ------------------------------------------------------------------------------
