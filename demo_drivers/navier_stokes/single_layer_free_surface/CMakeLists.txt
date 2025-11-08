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
message(VERBOSE "Entered single_layer_free_surface subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(single_layer_free_surface C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()

include(CTest)

oomph_add_executable(
  NAME elastic_single_layer
  SOURCES elastic_single_layer.cc
  LIBRARIES oomph::navier_stokes oomph::fluid_interface oomph::solid
            oomph::constitutive oomph::meshes oomph::generic)

oomph_add_executable(
  NAME spine_single_layer
  SOURCES spine_single_layer.cc
  LIBRARIES oomph::navier_stokes oomph::fluid_interface oomph::meshes
            oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME navier_stokes.single_layer_free_surface
  DEPENDS_ON elastic_single_layer spine_single_layer
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving single_layer_free_surface subdirectory")
# ------------------------------------------------------------------------------
