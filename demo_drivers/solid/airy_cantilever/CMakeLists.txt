# ------------------------------------------------------------------------------
# OOMPH auto-inserted: warn once if user configures inside an existing build.
if(EXISTS "${CMAKE_CURRENT_BINARY_DIR}/CMakeCache.txt"
   AND CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  message(
    WARNING
      "You're re-configuring an existing build tree. For a clean configuration, "
      "delete the whole build folder and start with a fresh directory.")
endif()
# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered airy_cantilever subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(airy_cantilever C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()

include(CTest)

oomph_add_executable(
  NAME airy_cantilever
  SOURCES airy_cantilever.cc
  LIBRARIES oomph::solid oomph::constitutive oomph::meshes oomph::generic)

oomph_add_executable(
  NAME airy_cantilever2_adapt
  SOURCES airy_cantilever2.cc
  LIBRARIES oomph::solid oomph::constitutive oomph::meshes oomph::generic
  CXX_DEFINITIONS REFINE)

oomph_add_executable(
  NAME airy_cantilever2_noadapt
  SOURCES airy_cantilever2.cc
  LIBRARIES oomph::solid oomph::constitutive oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME solid.airy_cantilever
  DEPENDS_ON airy_cantilever airy_cantilever2_adapt airy_cantilever2_noadapt
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving airy_cantilever subdirectory")
# ------------------------------------------------------------------------------
