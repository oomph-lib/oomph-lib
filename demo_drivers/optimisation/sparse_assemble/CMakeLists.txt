# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered sparse_assemble subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(sparse_assemble C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME sparse_assemble_test
  SOURCES sparse_assemble_test.cc
  LIBRARIES oomph::navier_stokes oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME optimisation.sparse_assemble
  DEPENDS_ON sparse_assemble_test
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata
  LABELS optimisation sparse_assemble)

message(VERBOSE "Leaving sparse_assemble subdirectory")
# ------------------------------------------------------------------------------
