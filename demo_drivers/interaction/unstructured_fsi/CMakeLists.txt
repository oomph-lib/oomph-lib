# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered unstructured_fsi subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(unstructured_fsi C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME unstructured_two_d_fsi
  SOURCES unstructured_two_d_fsi.cc
  LIBRARIES oomph::solid oomph::constitutive oomph::navier_stokes oomph::meshes
            oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME interaction.unstructured_fsi
  DEPENDS_ON unstructured_two_d_fsi
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh
             validata
             fluid.fig
             fluid.fig.1.ele
             fluid.fig.1.node
             fluid.fig.1.poly
             fluid.fig.poly
             solid.fig
             solid.fig.1.ele
             solid.fig.1.node
             solid.fig.1.poly
             solid.fig.poly)

message(VERBOSE "Leaving unstructured_fsi subdirectory")
# ------------------------------------------------------------------------------
