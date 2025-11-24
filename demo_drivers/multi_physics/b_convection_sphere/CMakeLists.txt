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
message(VERBOSE "Entered b_convection_sphere subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(b_convection_sphere C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()

include(CTest)

oomph_add_executable(
  NAME b_convection_sphere
  SOURCES half_rectangle_with_hole_mesh.h b_convection_sphere.cc
          axisym_buoyant_navier_stokes.h
  LIBRARIES oomph::axisym_advection_diffusion oomph::axisym_navier_stokes
            oomph::meshes oomph::generic)

oomph_add_executable(
  NAME multi_domain_b_convect_sph
  SOURCES multi_domain_b_convect_sph.cc half_rectangle_with_hole_mesh.h
          multi_domain_axisym_boussinesq_elements.h
  LIBRARIES oomph::axisym_advection_diffusion oomph::axisym_navier_stokes
            oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME multi_physics.b_convection_sphere
  DEPENDS_ON b_convection_sphere multi_domain_b_convect_sph
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving b_convection_sphere subdirectory")
# ------------------------------------------------------------------------------
