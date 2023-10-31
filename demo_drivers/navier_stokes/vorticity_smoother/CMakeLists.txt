# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered vorticity_smoother subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(vorticity_smoother C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME vorticity_smoother_validation
  SOURCES vorticity_smoother_validation.cc
  LIBRARIES oomph::navier_stokes oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME navier_stokes.vorticity_smoother
  DEPENDS_ON vorticity_smoother_validation
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata
  LABELS navier_stokes vorticity_smoother)

message(VERBOSE "Leaving vorticity_smoother subdirectory")
# ------------------------------------------------------------------------------
