# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered compressed_square subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(compressed_square C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME compressed_square
  SOURCES compressed_square.cc
  LIBRARIES oomph::solid oomph::constitutive oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME solid.compressed_square
  DEPENDS_ON compressed_square
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata
  LABELS solid compressed_square)

message(VERBOSE "Leaving compressed_square subdirectory")
# ------------------------------------------------------------------------------
