# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered vmtk_fsi subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(vmtk_fsi C CXX Fortran)
if(NOT OOMPH_BUILD_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME vmtk_fsi
  SOURCES vmtk_fsi.cc
  LIBRARIES oomph::solid oomph::constitutive oomph::navier_stokes oomph::meshes
            oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME interaction.vmtk_fsi
  TARGET_DEPENDENCIES vmtk_fsi
  EXTRA_REQUIRES
    fluid_iliac.1.ele
    fluid_iliac.1.face
    fluid_iliac.1.node
    fluid_iliac.poly
    iliac.xda
    solid_iliac.1.ele
    solid_iliac.1.face
    solid_iliac.1.node
    solid_iliac.poly
  LABELS interaction vmtk_fsi)

message(VERBOSE "Leaving vmtk_fsi subdirectory")
# ------------------------------------------------------------------------------
