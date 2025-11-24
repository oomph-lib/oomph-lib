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
message(VERBOSE "Entered sphere_scattering subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(sphere_scattering C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()

include(CTest)

oomph_add_executable(
  NAME sphere_scattering
  SOURCES sphere_scattering.cc
  LIBRARIES oomph::fourier_decomposed_helmholtz oomph::meshes oomph::generic)

oomph_add_executable(
  NAME unstructured_adaptive_sphere_scattering
  SOURCES unstructured_sphere_scattering.cc
  LIBRARIES oomph::fourier_decomposed_helmholtz oomph::meshes oomph::generic
  CXX_DEFINITIONS ADAPTIVE)

oomph_add_executable(
  NAME unstructured_sphere_scattering
  SOURCES unstructured_sphere_scattering.cc
  LIBRARIES oomph::fourier_decomposed_helmholtz oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME fourier_decomposed_helmholtz.sphere_scattering
  DEPENDS_ON sphere_scattering unstructured_adaptive_sphere_scattering
             unstructured_sphere_scattering
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving sphere_scattering subdirectory")
# ------------------------------------------------------------------------------
