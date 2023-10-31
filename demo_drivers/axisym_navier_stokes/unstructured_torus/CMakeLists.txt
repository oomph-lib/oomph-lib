# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered unstructured_torus subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(unstructured_torus C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME pressure_driven_torus
  SOURCES pressure_driven_torus.cc
  LIBRARIES oomph::solid oomph::constitutive oomph::axisym_navier_stokes
            oomph::meshes oomph::generic)

oomph_add_executable(
  NAME unstructured_torus
  SOURCES unstructured_torus.cc
  LIBRARIES oomph::axisym_navier_stokes oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME axisym_navier_stokes.unstructured_torus
  DEPENDS_ON pressure_driven_torus unstructured_torus
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata
  LABELS axisym_navier_stokes unstructured_torus)

message(VERBOSE "Leaving unstructured_torus subdirectory")
# ------------------------------------------------------------------------------
