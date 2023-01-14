# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered vmtk_fluid subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(vmtk_fluid C CXX Fortran)
if(NOT OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY)
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
  TARGET_DEPENDENCIES vmtk_fluid
  VALIDATE_SH_ARGS ${OOMPH_ROOT_DIR}
  EXTRA_REQUIRES fluid_iliac.1.ele fluid_iliac.1.face fluid_iliac.1.node
                 fluid_iliac.poly iliac.xda
  LABELS navier_stokes vmtk_fluid)

message(VERBOSE "Leaving vmtk_fluid subdirectory")
# ------------------------------------------------------------------------------
