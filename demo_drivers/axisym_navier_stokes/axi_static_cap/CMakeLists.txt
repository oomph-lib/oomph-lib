# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered axi_static_cap subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(axi_static_cap C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME axi_static_cap
  SOURCES axi_static_cap.cc
  LIBRARIES oomph::axisym_navier_stokes oomph::fluid_interface oomph::solid
            oomph::constitutive oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME axisym_navier_stokes.axi_static_cap
  DEPENDS_ON axi_static_cap
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata
  LABELS axisym_navier_stokes axi_static_cap)

message(VERBOSE "Leaving axi_static_cap subdirectory")
# ------------------------------------------------------------------------------
