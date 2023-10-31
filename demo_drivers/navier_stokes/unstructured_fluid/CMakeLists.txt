# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered unstructured_fluid subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(unstructured_fluid C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME unstructured_two_d_fluid
  SOURCES unstructured_two_d_fluid.cc
  LIBRARIES oomph::solid oomph::constitutive oomph::navier_stokes oomph::meshes
            oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME navier_stokes.unstructured_fluid
  DEPENDS_ON unstructured_two_d_fluid
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh
             validata
             fluid.fig
             fluid.fig.1.ele
             fluid.fig.1.node
             fluid.fig.1.poly
             fluid.fig.poly
  LABELS navier_stokes unstructured_fluid)

message(VERBOSE "Leaving unstructured_fluid subdirectory")
# ------------------------------------------------------------------------------
