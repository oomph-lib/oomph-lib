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
message(VERBOSE "Entered jeffery_hamel subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(jeffery_hamel C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()

include(CTest)

oomph_add_executable(
  NAME jeffery_hamel
  SOURCES jeffery_hamel.cc
          refineable_r_mesh.h
          refineable_streamfunction_elements.h
          streamfunction_elements.h
          streamfunction_include.h
          jh_includes.h
          polar_streamfunction_traction_elements.h
  LIBRARIES oomph::polar_navier_stokes oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME polar_navier_stokes.jeffery_hamel
  DEPENDS_ON jeffery_hamel
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving jeffery_hamel subdirectory")
# ------------------------------------------------------------------------------
