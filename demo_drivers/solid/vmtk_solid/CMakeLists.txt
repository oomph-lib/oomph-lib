# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered vmtk_solid subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(vmtk_solid C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME vmtk_solid
  SOURCES vmtk_solid.cc
  LIBRARIES oomph::solid oomph::constitutive oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME solid.vmtk_solid
  DEPENDS_ON vmtk_solid
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata solid_iliac.1.ele solid_iliac.1.face
             solid_iliac.1.node solid_iliac.poly)

message(VERBOSE "Leaving vmtk_solid subdirectory")
# ------------------------------------------------------------------------------
