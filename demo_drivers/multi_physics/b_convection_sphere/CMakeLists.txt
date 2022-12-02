# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered b_convection_sphere subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(b_convection_sphere C CXX Fortran)
if(NOT BUILD_DEMO_DRIVERS_WITH_LIBRARY)
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
  TARGET_DEPENDENCIES b_convection_sphere multi_domain_b_convect_sph
  EXTRA_REQUIRES
  LABELS multi_physics b_convection_sphere)

message(VERBOSE "Leaving b_convection_sphere subdirectory")
# ------------------------------------------------------------------------------
