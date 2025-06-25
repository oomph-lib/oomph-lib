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
message(VERBOSE "Entered axisym_heat_sphere subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(axisym_heat_sphere C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()

include(CTest)

oomph_add_executable(
  NAME axisym_heat_sphere
  SOURCES half_rectangle_with_hole_mesh.h axisym_heat_sphere.cc
          axisym_buoyant_navier_stokes.h
  LIBRARIES oomph::axisym_advection_diffusion oomph::axisym_navier_stokes
            oomph::meshes oomph::generic)

oomph_add_executable(
  NAME multi_domain_axisym_heat_sphere
  SOURCES multi_domain_axisym_heat_sphere.cc half_rectangle_with_hole_mesh.h
          multi_domain_axisym_boussinesq_elements.h
  LIBRARIES oomph::axisym_advection_diffusion oomph::axisym_navier_stokes
            oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME multi_physics.axisym_heat_sphere
  DEPENDS_ON axisym_heat_sphere multi_domain_axisym_heat_sphere
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving axisym_heat_sphere subdirectory")
# ------------------------------------------------------------------------------
