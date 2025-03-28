# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered static_cap subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(static_cap C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME static_single_layer
  SOURCES static_single_layer.cc
  LIBRARIES oomph::navier_stokes oomph::fluid_interface oomph::solid
            oomph::constitutive oomph::meshes oomph::generic)

oomph_add_executable(
  NAME static_two_layer
  SOURCES static_two_layer.cc
  LIBRARIES oomph::navier_stokes oomph::fluid_interface oomph::solid
            oomph::constitutive oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME navier_stokes.static_cap
  DEPENDS_ON static_single_layer static_two_layer
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving static_cap subdirectory")
# ------------------------------------------------------------------------------
