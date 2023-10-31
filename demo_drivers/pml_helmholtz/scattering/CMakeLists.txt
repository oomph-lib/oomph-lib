# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered scattering subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(scattering C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME unstructured_two_d_helmholtz
  SOURCES unstructured_two_d_helmholtz.cc
  LIBRARIES oomph::pml_helmholtz oomph::meshes oomph::generic)

oomph_add_executable(
  NAME unstructured_two_d_helmholtz_adapt
  SOURCES unstructured_two_d_helmholtz.cc
  LIBRARIES oomph::pml_helmholtz oomph::meshes oomph::generic
  CXX_DEFINITIONS ADAPTIVE)

oomph_add_executable(
  NAME unstructured_two_d_helmholtz_scattering
  SOURCES unstructured_two_d_helmholtz_scattering.cc
  LIBRARIES oomph::pml_helmholtz oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME pml_helmholtz.scattering
  DEPENDS_ON unstructured_two_d_helmholtz unstructured_two_d_helmholtz_adapt
             unstructured_two_d_helmholtz_scattering
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata
  LABELS pml_helmholtz scattering)

message(VERBOSE "Leaving scattering subdirectory")
# ------------------------------------------------------------------------------
