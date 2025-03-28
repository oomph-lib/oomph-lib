# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered jeffery_hamel subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(jeffery_hamel C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME jeffery_hamel
  SOURCES jeffery_hamel.cc
          refineable_r_mesh.h
          refineable_streamfunction_elements.h
          streamfunction_elements.h
          streamfunction_include.h
          jh_includes.h
          polar_streamfunction_traction_elements.h
  LIBRARIES oomph::polar_navier_stokes oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME polar_navier_stokes.jeffery_hamel
  DEPENDS_ON jeffery_hamel
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving jeffery_hamel subdirectory")
# ------------------------------------------------------------------------------
