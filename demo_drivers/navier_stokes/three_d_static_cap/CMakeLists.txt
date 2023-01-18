# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered three_d_static_cap subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(three_d_static_cap C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME 3d_static_cap
  SOURCES 3d_static_cap.cc
  LIBRARIES oomph::navier_stokes oomph::fluid_interface oomph::meshes
            oomph::generic)

oomph_add_executable(
  NAME 3d_static_cap_elastic
  SOURCES 3d_static_cap_elastic.cc
  LIBRARIES oomph::navier_stokes oomph::fluid_interface oomph::solid
            oomph::constitutive oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME navier_stokes.three_d_static_cap
  DEPENDS_ON 3d_static_cap 3d_static_cap_elastic
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata
  LABELS navier_stokes three_d_static_cap)

message(VERBOSE "Leaving three_d_static_cap subdirectory")
# ------------------------------------------------------------------------------
