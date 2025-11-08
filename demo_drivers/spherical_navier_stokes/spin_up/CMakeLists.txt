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
message(VERBOSE "Entered spin_up subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(spin_up C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()

include(CTest)

oomph_add_executable(
  NAME refineable_spin_up
  SOURCES refineable_spin_up.cc
  LIBRARIES oomph::spherical_navier_stokes oomph::meshes oomph::generic)

oomph_add_executable(
  NAME refineable_spin_up_cyl
  SOURCES refineable_spin_up_cyl.cc
  LIBRARIES oomph::axisym_navier_stokes oomph::meshes oomph::generic)

oomph_add_executable(
  NAME spin_up
  SOURCES spin_up.cc
  LIBRARIES oomph::spherical_navier_stokes oomph::meshes oomph::generic)

oomph_add_executable(
  NAME spin_up_cyl
  SOURCES spin_up_cyl.cc
  LIBRARIES oomph::axisym_navier_stokes oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME spherical_navier_stokes.spin_up
  DEPENDS_ON refineable_spin_up refineable_spin_up_cyl spin_up spin_up_cyl
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving spin_up subdirectory")
# ------------------------------------------------------------------------------
