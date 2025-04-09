# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered basic_ode subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(basic_ode C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME basic_ode
  SOURCES basic_ode.cc
  LIBRARIES oomph::ode oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME ode.basic_ode
  DEPENDS_ON basic_ode
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh zeros
  SILENCE_MISSING_VALIDATA_WARNING)

message(VERBOSE "Leaving basic_ode subdirectory")
# ------------------------------------------------------------------------------
