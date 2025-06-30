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
message(VERBOSE "Entered two_d_linear_wave subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(two_d_linear_wave C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()

include(CTest)

oomph_add_executable(
  NAME two_d_linear_wave
  SOURCES two_d_linear_wave.cc
  LIBRARIES oomph::linear_wave oomph::meshes oomph::generic)

oomph_add_executable(
  NAME two_d_linear_wave_flux
  SOURCES two_d_linear_wave_flux.cc
  LIBRARIES oomph::linear_wave oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME linear_wave.two_d_linear_wave
  DEPENDS_ON two_d_linear_wave two_d_linear_wave_flux
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving two_d_linear_wave subdirectory")
# ------------------------------------------------------------------------------
