# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered elastic_annulus subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(elastic_annulus C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME adaptive_time_harmonic_elastic_annulus
  SOURCES time_harmonic_elastic_annulus.cc
  LIBRARIES oomph::time_harmonic_linear_elasticity oomph::meshes oomph::generic
  CXX_DEFINITIONS ADAPTIVE)

oomph_add_executable(
  NAME adaptive_unstructured_time_harmonic_elastic_annulus
  SOURCES unstructured_time_harmonic_elastic_annulus.cc
  LIBRARIES oomph::time_harmonic_linear_elasticity oomph::meshes oomph::generic
  CXX_DEFINITIONS ADAPTIVE)

oomph_add_executable(
  NAME time_harmonic_elastic_annulus
  SOURCES time_harmonic_elastic_annulus.cc
  LIBRARIES oomph::time_harmonic_linear_elasticity oomph::meshes oomph::generic)

oomph_add_executable(
  NAME unstructured_time_harmonic_elastic_annulus
  SOURCES unstructured_time_harmonic_elastic_annulus.cc
  LIBRARIES oomph::time_harmonic_linear_elasticity oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME time_harmonic_linear_elasticity.elastic_annulus
  DEPENDS_ON adaptive_time_harmonic_elastic_annulus
             adaptive_unstructured_time_harmonic_elastic_annulus
             time_harmonic_elastic_annulus
             unstructured_time_harmonic_elastic_annulus
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving elastic_annulus subdirectory")
# ------------------------------------------------------------------------------
