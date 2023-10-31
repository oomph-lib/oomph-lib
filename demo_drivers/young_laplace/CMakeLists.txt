# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered young_laplace subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(young_laplace C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME barrel
  SOURCES barrel.cc
  LIBRARIES oomph::young_laplace oomph::meshes oomph::generic)

oomph_add_executable(
  NAME refineable_t_junction
  SOURCES refineable_t_junction.cc
  LIBRARIES oomph::young_laplace oomph::meshes oomph::generic
  CXX_OPTIONS -g)

oomph_add_executable(
  NAME refineable_young_laplace
  SOURCES refineable_young_laplace.cc common_young_laplace_stuff.h
  LIBRARIES oomph::young_laplace oomph::meshes oomph::generic)

oomph_add_executable(
  NAME spherical_cap_in_cylinder
  SOURCES spherical_cap_in_cylinder.cc common_young_laplace_stuff.h
  LIBRARIES oomph::young_laplace oomph::meshes oomph::generic)

oomph_add_executable(
  NAME young_laplace
  SOURCES young_laplace.cc common_young_laplace_stuff.h
  LIBRARIES oomph::young_laplace oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME young_laplace.young_laplace
  DEPENDS_ON barrel refineable_t_junction refineable_young_laplace
             spherical_cap_in_cylinder young_laplace
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata
  LABELS young_laplace)

message(VERBOSE "Leaving young_laplace subdirectory")
# ------------------------------------------------------------------------------
