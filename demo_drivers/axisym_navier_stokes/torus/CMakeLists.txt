# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered torus subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(torus C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME torus
  SOURCES torus.cc
  LIBRARIES oomph::navier_stokes oomph::axisym_navier_stokes oomph::meshes
            oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME axisym_navier_stokes.torus
  DEPENDS_ON torus
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata
  LABELS axisym_navier_stokes torus)

message(VERBOSE "Leaving torus subdirectory")
# ------------------------------------------------------------------------------
