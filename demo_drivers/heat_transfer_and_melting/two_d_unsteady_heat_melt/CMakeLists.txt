# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered two_d_unsteady_heat_melt subdirectory")

cmake_minimum_required(VERSION 3.16..3.20 FATAL_ERROR)
project(two_d_unsteady_heat_melt C CXX Fortran)
if(NOT BUILDING_DEMO_DRIVERS_WITH_LIBRARY)
  find_package(oomphlib REQUIRED)
endif()
include(CTest)

oomph_add_executable(
  NAME circular_boulder_melt
  SOURCES circular_boulder_melt.cc
  LIBRARIES oomph::unsteady_heat oomph::constitutive oomph::solid
  CXX_STANDARD 14)

oomph_add_executable(
  NAME circular_boulder_solar_radiation
  SOURCES circular_boulder_solar_radiation.cc
  LIBRARIES oomph::unsteady_heat oomph::constitutive oomph::solid
  CXX_STANDARD 14)

oomph_add_executable(
  NAME heated_linear_solid_contact_with_gravity
  SOURCES heated_linear_solid_contact_with_gravity.cc contact_elements.h
  LIBRARIES oomph::unsteady_heat oomph::linear_elasticity oomph::constitutive
            oomph::solid
  CXX_STANDARD 14)

oomph_add_executable(
  NAME linear_solid_contact_with_gravity_structured
  SOURCES linear_solid_contact_with_gravity.cc
  LIBRARIES oomph::linear_elasticity oomph::constitutive oomph::solid
  CXX_DEFINITIONS STRUCTURED_MESH
  CXX_STANDARD 14)

oomph_add_executable(
  NAME linear_solid_contact_with_gravity_unstructured
  SOURCES linear_solid_contact_with_gravity.cc
  LIBRARIES oomph::linear_elasticity oomph::constitutive oomph::solid
  CXX_STANDARD 14)

oomph_add_executable(
  NAME melt
  SOURCES melt.cc heat_transfer_and_melt_elements.h
  LIBRARIES oomph::constitutive oomph::solid oomph::unsteady_heat
  CXX_STANDARD 14)

oomph_add_executable(
  NAME pretend_melt
  SOURCES pretend_melt.cc unsteady_heat_flux_pseudo_melt_elements.h
  LIBRARIES oomph::unsteady_heat
  CXX_STANDARD 14)

oomph_add_executable(
  NAME solid_contact
  SOURCES solid_contact.cc
  LIBRARIES oomph::solid oomph::constitutive
  CXX_STANDARD 14)

oomph_add_executable(
  NAME solid_contact_with_gravity_structured
  SOURCES solid_contact_with_gravity.cc contact_elements.h
  LIBRARIES oomph::constitutive oomph::solid
  CXX_DEFINITIONS STRUCTURED_MESH
  CXX_STANDARD 14)

oomph_add_executable(
  NAME solid_contact_with_gravity_unstructured
  SOURCES solid_contact_with_gravity.cc contact_elements.h
  LIBRARIES oomph::constitutive oomph::solid
  CXX_STANDARD 14)

oomph_add_executable(
  NAME spring_contact
  SOURCES spring_contact.cc
  LIBRARIES oomph::generic
  CXX_STANDARD 14)

oomph_add_executable(
  NAME spring_contact2
  SOURCES spring_contact2.cc
  LIBRARIES oomph::generic
  CXX_STANDARD 14)

oomph_add_executable(
  NAME stefan_boltzmann
  SOURCES stefan_boltzmann.cc temporary_stefan_boltzmann_elements.h
  LIBRARIES oomph::constitutive oomph::solid oomph::unsteady_heat
  CXX_STANDARD 14)

oomph_add_executable(
  NAME stefan_boltzmann_melt
  SOURCES stefan_boltzmann_melt.cc heat_transfer_and_melt_elements.h
          temporary_stefan_boltzmann_elements.h
  LIBRARIES oomph::constitutive oomph::solid oomph::unsteady_heat
  CXX_STANDARD 14)

oomph_add_executable(
  NAME two_d_unsteady_heat
  SOURCES two_d_unsteady_heat.cc
  LIBRARIES oomph::unsteady_heat
  CXX_STANDARD 14)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME heat_transfer_and_melting.two_d_unsteady_heat_melt
  TARGET_DEPENDENCIES
    circular_boulder_melt
    circular_boulder_solar_radiation
    heated_linear_solid_contact_with_gravity
    linear_solid_contact_with_gravity_structured
    linear_solid_contact_with_gravity_unstructured
    melt
    pretend_melt
    solid_contact
    solid_contact_with_gravity_structured
    solid_contact_with_gravity_unstructured
    spring_contact
    spring_contact2
    stefan_boltzmann
    stefan_boltzmann_melt
    two_d_unsteady_heat
  EXTRA_REQUIRES
  LABELS "heat_transfer_and_melting;two_d_unsteady_heat_melt")

# This test contains a huge number of subtests and as such needs a larger amount
# of time than the default to complete
set_tests_properties(heat_transfer_and_melting.two_d_unsteady_heat_melt
                     PROPERTIES TIMEOUT 10000)

message(VERBOSE "Leaving two_d_unsteady_heat_melt subdirectory")
# ------------------------------------------------------------------------------
