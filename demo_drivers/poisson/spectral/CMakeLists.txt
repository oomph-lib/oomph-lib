# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered spectral subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(spectral C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME one_d_spectral
  SOURCES one_d_spectral.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

oomph_add_executable(
  NAME three_d_spectral_adapt
  SOURCES three_d_spectral_adapt.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

oomph_add_executable(
  NAME two_d_spectral_adapt
  SOURCES two_d_spectral_adapt.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME poisson.spectral
  DEPENDS_ON one_d_spectral three_d_spectral_adapt two_d_spectral_adapt
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata
  LABELS poisson spectral)

message(VERBOSE "Leaving spectral subdirectory")
# ------------------------------------------------------------------------------
