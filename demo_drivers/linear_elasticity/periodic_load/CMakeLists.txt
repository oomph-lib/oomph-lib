# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered periodic_load subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(periodic_load C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME periodic_load
  SOURCES periodic_load.cc
  LIBRARIES oomph::linear_elasticity oomph::meshes oomph::generic)

oomph_add_executable(
  NAME refineable_periodic_load
  SOURCES refineable_periodic_load.cc
  LIBRARIES oomph::linear_elasticity oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME linear_elasticity.periodic_load
  DEPENDS_ON periodic_load refineable_periodic_load
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving periodic_load subdirectory")
# ------------------------------------------------------------------------------
