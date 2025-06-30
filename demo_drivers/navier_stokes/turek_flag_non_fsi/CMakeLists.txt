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
message(VERBOSE "Entered turek_flag_non_fsi subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(turek_flag_non_fsi C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()

include(CTest)

oomph_add_executable(
  NAME turek_flag_non_fsi_alg
  SOURCES turek_flag_non_fsi.cc
  LIBRARIES oomph::navier_stokes oomph::meshes oomph::generic)

oomph_add_executable(
  NAME turek_flag_non_fsi_macro
  SOURCES turek_flag_non_fsi.cc
  LIBRARIES oomph::navier_stokes oomph::meshes oomph::generic
  CXX_DEFINITIONS USE_MACRO_ELEMENTS)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME navier_stokes.turek_flag_non_fsi
  DEPENDS_ON turek_flag_non_fsi_alg turek_flag_non_fsi_macro
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving turek_flag_non_fsi subdirectory")
# ------------------------------------------------------------------------------
