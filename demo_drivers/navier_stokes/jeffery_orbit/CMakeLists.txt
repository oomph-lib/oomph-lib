# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered jeffery_orbit subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(jeffery_orbit C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME jeffery_orbit
  SOURCES jeffery_orbit.cc my_taylor_hood_elements.h
  LIBRARIES oomph::solid oomph::constitutive oomph::navier_stokes
            oomph::rigid_body oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME navier_stokes.jeffery_orbit
  DEPENDS_ON jeffery_orbit
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata
  LABELS navier_stokes jeffery_orbit)

message(VERBOSE "Leaving jeffery_orbit subdirectory")
# ------------------------------------------------------------------------------
