# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered unstructured_adaptive_fsi subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(unstructured_adaptive_fsi C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME unstructured_adaptive_2d_fsi
  SOURCES unstructured_adaptive_2d_fsi.cc
  LIBRARIES oomph::solid oomph::constitutive oomph::navier_stokes oomph::meshes
            oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME interaction.unstructured_adaptive_fsi
  DEPENDS_ON unstructured_adaptive_2d_fsi
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving unstructured_adaptive_fsi subdirectory")
# ------------------------------------------------------------------------------
