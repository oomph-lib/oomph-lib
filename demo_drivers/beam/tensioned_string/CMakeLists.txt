# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered tensioned_string subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(tensioned_string C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME tensioned_string
  SOURCES tensioned_string.cc
  LIBRARIES oomph::beam oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME beam.tensioned_string
  DEPENDS_ON tensioned_string
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata
  LABELS beam tensioned_string)

message(VERBOSE "Leaving tensioned_string subdirectory")
# ------------------------------------------------------------------------------
