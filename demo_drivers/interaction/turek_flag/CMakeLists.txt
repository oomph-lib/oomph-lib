# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered turek_flag subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(turek_flag C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME turek_flag
  SOURCES turek_flag.cc
  LIBRARIES oomph::solid oomph::constitutive oomph::navier_stokes oomph::meshes
            oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME interaction.turek_flag
  DEPENDS_ON turek_flag
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving turek_flag subdirectory")
# ------------------------------------------------------------------------------
