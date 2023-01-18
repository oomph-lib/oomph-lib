# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered elastic_bretherton subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(elastic_bretherton C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME elastic_bretherton
  SOURCES elastic_bretherton.cc
  LIBRARIES oomph::beam oomph::navier_stokes oomph::fluid_interface
            oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME interaction.elastic_bretherton
  DEPENDS_ON elastic_bretherton
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata
  LABELS interaction elastic_bretherton)

message(VERBOSE "Leaving elastic_bretherton subdirectory")
# ------------------------------------------------------------------------------
