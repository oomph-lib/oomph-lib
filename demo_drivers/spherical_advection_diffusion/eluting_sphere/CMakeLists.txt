# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered eluting_sphere subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(eluting_sphere C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME eluting_sphere
  SOURCES eluting_sphere.cc
  LIBRARIES oomph::spherical_advection_diffusion oomph::meshes oomph::generic)

oomph_add_executable(
  NAME eluting_sphere_adapt
  SOURCES eluting_sphere_adapt.cc
  LIBRARIES oomph::spherical_advection_diffusion oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME spherical_advection_diffusion.eluting_sphere
  DEPENDS_ON eluting_sphere eluting_sphere_adapt
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving eluting_sphere subdirectory")
# ------------------------------------------------------------------------------
