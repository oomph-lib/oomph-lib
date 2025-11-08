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
message(VERBOSE "Entered jeffery_orbit subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(jeffery_orbit C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()

include(CTest)

oomph_add_executable(
  NAME jeffery_orbit
  SOURCES jeffery_orbit.cc my_taylor_hood_elements.h
  LIBRARIES oomph::solid oomph::constitutive oomph::navier_stokes
            oomph::rigid_body oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME navier_stokes.jeffery_orbit
  DEPENDS_ON jeffery_orbit
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving jeffery_orbit subdirectory")
# ------------------------------------------------------------------------------
