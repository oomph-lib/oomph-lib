# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered navier_stokes subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(navier_stokes C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME circular_driven_cavity
  SOURCES circular_driven_cavity_driver.cc validate.sh
  LIBRARIES oomph::navier_stokes oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME self_test.navier_stokes
  DEPENDS_ON circular_driven_cavity
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving navier_stokes subdirectory")
# ------------------------------------------------------------------------------
