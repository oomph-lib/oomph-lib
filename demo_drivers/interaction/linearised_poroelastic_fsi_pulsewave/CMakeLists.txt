# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered linearised_poroelastic_fsi_pulsewave subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(linearised_poroelastic_fsi_pulsewave C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME linearised_poroelastic_fsi_pulsewave
  SOURCES linearised_poroelastic_fsi_pulsewave.cc
  LIBRARIES oomph::axisym_navier_stokes oomph::axisym_poroelasticity
            oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME interaction.linearised_poroelastic_fsi_pulsewave
  DEPENDS_ON linearised_poroelastic_fsi_pulsewave
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata
  LABELS interaction linearised_poroelastic_fsi_pulsewave)

message(VERBOSE "Leaving linearised_poroelastic_fsi_pulsewave subdirectory")
# ------------------------------------------------------------------------------
