# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered elastic_poisson subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(elastic_poisson C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME elastic_poisson
  SOURCES elastic_poisson.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME poisson.elastic_poisson
  DEPENDS_ON elastic_poisson
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving elastic_poisson subdirectory")
# ------------------------------------------------------------------------------
