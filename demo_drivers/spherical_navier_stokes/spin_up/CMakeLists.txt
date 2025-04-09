# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered spin_up subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(spin_up C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME refineable_spin_up
  SOURCES refineable_spin_up.cc
  LIBRARIES oomph::spherical_navier_stokes oomph::meshes oomph::generic)

oomph_add_executable(
  NAME refineable_spin_up_cyl
  SOURCES refineable_spin_up_cyl.cc
  LIBRARIES oomph::axisym_navier_stokes oomph::meshes oomph::generic)

oomph_add_executable(
  NAME spin_up
  SOURCES spin_up.cc
  LIBRARIES oomph::spherical_navier_stokes oomph::meshes oomph::generic)

oomph_add_executable(
  NAME spin_up_cyl
  SOURCES spin_up_cyl.cc
  LIBRARIES oomph::axisym_navier_stokes oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME spherical_navier_stokes.spin_up
  DEPENDS_ON refineable_spin_up refineable_spin_up_cyl spin_up spin_up_cyl
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving spin_up subdirectory")
# ------------------------------------------------------------------------------
