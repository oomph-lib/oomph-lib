# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered two_d_variable_diff subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(two_d_variable_diff C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME two_d_variable_diff_adapt
  SOURCES two_d_variable_diff_adapt.cc
  LIBRARIES oomph::advection_diffusion oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME advection_diffusion.two_d_variable_diff
  DEPENDS_ON two_d_variable_diff_adapt
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving two_d_variable_diff subdirectory")
# ------------------------------------------------------------------------------
