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
message(VERBOSE "Entered bretherton subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(bretherton C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()

include(CTest)

oomph_add_executable(
  NAME bretherton
  SOURCES bretherton.cc
  LIBRARIES oomph::navier_stokes oomph::fluid_interface oomph::meshes
            oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME navier_stokes.bretherton
  DEPENDS_ON bretherton
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

# This sometimes (rather bizarrely) just times out on macOS. It shouldn't take
# more than 10 seconds to run though...
if(APPLE)
  set_tests_properties(navier_stokes.bretherton PROPERTIES DISABLED TRUE)
endif()

message(VERBOSE "Leaving bretherton subdirectory")
# ------------------------------------------------------------------------------
