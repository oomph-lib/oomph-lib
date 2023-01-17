# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered disk_compression subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(disk_compression C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME disk_compression
  SOURCES disk_compression.cc
  LIBRARIES oomph::solid oomph::constitutive oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME solid.disk_compression
  TARGET_DEPENDENCIES disk_compression
  VALIDATE_SH_ARGS ${OOMPH_ROOT_DIR}
  EXTRA_REQUIRES
  LABELS solid disk_compression)

message(VERBOSE "Leaving disk_compression subdirectory")
# ------------------------------------------------------------------------------
