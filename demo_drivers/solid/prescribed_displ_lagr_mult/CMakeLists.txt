# ------------------------------------------------------------------------------
# OOMPH auto-inserted: warn once if user configures inside an existing build.
if(EXISTS "${CMAKE_CURRENT_BINARY_DIR}/CMakeCache.txt"
   AND CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  message(
    WARNING
      "You're re-configuring an existing build tree. For a clean configuration, "
      "delete the whole build folder and start with a fresh directory.")
endif()
# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered prescribed_displ_lagr_mult subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(prescribed_displ_lagr_mult C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()

include(CTest)

oomph_add_executable(
  NAME prescribed_displ_lagr_mult
  SOURCES prescribed_displ_lagr_mult.cc
  LIBRARIES oomph::solid oomph::constitutive oomph::meshes oomph::generic)

oomph_add_executable(
  NAME prescribed_displ_lagr_mult2
  SOURCES prescribed_displ_lagr_mult2.cc
  LIBRARIES oomph::solid oomph::constitutive oomph::meshes oomph::generic)

oomph_add_executable(
  NAME prescribed_displ_lagr_mult_precond
  SOURCES prescribed_displ_lagr_mult_precond.cc
  LIBRARIES oomph::multi_physics oomph::solid oomph::constitutive
            oomph::navier_stokes oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME solid.prescribed_displ_lagr_mult
  DEPENDS_ON prescribed_displ_lagr_mult prescribed_displ_lagr_mult2
             prescribed_displ_lagr_mult_precond
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving prescribed_displ_lagr_mult subdirectory")
# ------------------------------------------------------------------------------
