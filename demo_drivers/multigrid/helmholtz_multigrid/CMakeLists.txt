# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered helmholtz_multigrid subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(helmholtz_multigrid C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME structured_cubic_point_source
  SOURCES structured_cubic_point_source.cc
  LIBRARIES oomph::pml_helmholtz oomph::meshes oomph::generic)

oomph_add_executable(
  NAME unstructured_two_d_helmholtz
  SOURCES unstructured_two_d_helmholtz.cc
  LIBRARIES oomph::pml_helmholtz oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME multigrid.helmholtz_multigrid
  DEPENDS_ON structured_cubic_point_source unstructured_two_d_helmholtz
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving helmholtz_multigrid subdirectory")
# ------------------------------------------------------------------------------
