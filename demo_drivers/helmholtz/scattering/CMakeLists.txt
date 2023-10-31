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
  NAME adaptive_scattering
  SOURCES scattering.cc
  LIBRARIES oomph::helmholtz oomph::meshes oomph::generic
  CXX_DEFINITIONS ADAPTIVE)

oomph_add_executable(
  NAME adaptive_unstructured_scattering
  SOURCES unstructured_scattering.cc
  LIBRARIES oomph::helmholtz oomph::meshes oomph::generic
  CXX_DEFINITIONS ADAPTIVE)

oomph_add_executable(
  NAME scattering
  SOURCES scattering.cc
  LIBRARIES oomph::helmholtz oomph::meshes oomph::generic)

oomph_add_executable(
  NAME unstructured_scattering
  SOURCES unstructured_scattering.cc
  LIBRARIES oomph::helmholtz oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME helmholtz.scattering
  DEPENDS_ON adaptive_scattering adaptive_unstructured_scattering scattering
             unstructured_scattering
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata
  LABELS helmholtz scattering)

message(VERBOSE "Leaving scattering subdirectory")
# ------------------------------------------------------------------------------
