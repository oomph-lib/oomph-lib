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
message(VERBOSE
        "Entered two_layer_interface_nonaxisym_perturbations subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(two_layer_interface_nonaxisym_perturbations C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()

include(CTest)

oomph_add_executable(
  NAME two_layer_interface_nonaxisym_perturbations
  SOURCES two_layer_interface_nonaxisym_perturbations.cc
  LIBRARIES oomph::axisym_navier_stokes oomph::navier_stokes
            oomph::fluid_interface oomph::oomph_crbond_bessel oomph::meshes
            oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME
    linearised_free_surface_axisym_navier_stokes.two_layer_interface_nonaxisym_perturbations
  DEPENDS_ON two_layer_interface_nonaxisym_perturbations
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE
        "Leaving two_layer_interface_nonaxisym_perturbations subdirectory")
# ------------------------------------------------------------------------------
