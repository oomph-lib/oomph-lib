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
message(VERBOSE "Entered simple_shear subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(simple_shear C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()

include(CTest)

oomph_add_executable(
  NAME simple_shear
  SOURCES simple_shear.cc precompiled_mesh.cc
  LIBRARIES oomph::solid oomph::constitutive oomph::generic oomph::meshes)

oomph_add_executable(
  NAME refineable_simple_shear
  SOURCES refineable_simple_shear.cc precompiled_mesh.cc
  LIBRARIES oomph::solid oomph::constitutive oomph::generic oomph::meshes)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME solid.simple_shear
  DEPENDS_ON simple_shear refineable_simple_shear
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

# This test needs more than the default allotted time to complete
set_tests_properties(solid.simple_shear PROPERTIES TIMEOUT 3600)

message(VERBOSE "Leaving simple_shear subdirectory")
# ------------------------------------------------------------------------------
