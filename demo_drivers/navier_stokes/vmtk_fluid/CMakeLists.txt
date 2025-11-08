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
message(VERBOSE "Entered vmtk_fluid subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(vmtk_fluid C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()

include(CTest)

oomph_add_executable(
  NAME vmtk_fluid
  SOURCES vmtk_fluid.cc
  LIBRARIES oomph::navier_stokes oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME navier_stokes.vmtk_fluid
  DEPENDS_ON vmtk_fluid
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh
             validata
             fluid_iliac.1.ele
             fluid_iliac.1.face
             fluid_iliac.1.node
             fluid_iliac.poly
             iliac.xda)

message(VERBOSE "Leaving vmtk_fluid subdirectory")
# ------------------------------------------------------------------------------
