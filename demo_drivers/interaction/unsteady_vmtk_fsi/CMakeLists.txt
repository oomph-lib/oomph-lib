# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered unsteady_vmtk_fsi subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(unsteady_vmtk_fsi C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME unsteady_vmtk_fsi
  SOURCES unsteady_vmtk_fsi.cc
  LIBRARIES oomph::solid
            oomph::constitutive
            oomph::navier_stokes
            oomph::poisson
            oomph::linear_elasticity
            oomph::meshes
            oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME interaction.unsteady_vmtk_fsi
  DEPENDS_ON unsteady_vmtk_fsi
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh
             validata
             fluid.1.ele
             fluid.1.face
             fluid.1.node
             fluid_quadratic_nodes.dat
             solid.1.ele
             solid.1.face
             solid.1.node
             solid_quadratic_nodes.dat
             boundary_enumeration.dat)

message(VERBOSE "Leaving unsteady_vmtk_fsi subdirectory")
# ------------------------------------------------------------------------------
