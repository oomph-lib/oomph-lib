# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered fibre subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(fibre C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME fibre
  SOURCES fibre.cc
  LIBRARIES oomph::axisym_navier_stokes oomph::fluid_interface oomph::meshes
            oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME axisym_navier_stokes.fibre
  DEPENDS_ON fibre
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving fibre subdirectory")
# ------------------------------------------------------------------------------
