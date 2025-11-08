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
message(VERBOSE "Entered free_boundary_poisson subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(free_boundary_poisson C CXX Fortran)

if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()

include(CTest)

oomph_add_executable(
  NAME algebraic_free_boundary_poisson
  SOURCES algebraic_free_boundary_poisson.cc circle_as_generalised_element.h
  LIBRARIES oomph::poisson oomph::constitutive oomph::meshes oomph::generic)

oomph_add_executable(
  NAME circle
  SOURCES circle.cc circle.h
  LIBRARIES oomph::meshes oomph::generic)

oomph_add_executable(
  NAME doc_sparse_macro_node_update
  SOURCES doc_sparse_macro_node_update.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

oomph_add_executable(
  NAME elastic_poisson
  SOURCES elastic_poisson.cc circle_as_generalised_element.h
  LIBRARIES oomph::poisson oomph::solid oomph::constitutive oomph::meshes
            oomph::generic)

oomph_add_executable(
  NAME free_boundary_poisson_with_elastic_mesh_update
  SOURCES elastic_mesh_update.cc
  LIBRARIES oomph::poisson oomph::solid oomph::constitutive oomph::meshes
            oomph::generic)

oomph_add_executable(
  NAME geom_object_element
  SOURCES geom_object_element.cc circle_as_generalised_element.h
  LIBRARIES oomph::meshes oomph::generic)

oomph_add_executable(
  NAME macro_element_free_boundary_poisson
  SOURCES macro_element_free_boundary_poisson.cc circle_as_generalised_element.h
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

oomph_add_executable(
  NAME macro_element_free_boundary_poisson_non_ref
  SOURCES macro_element_free_boundary_poisson_non_ref.cc
          circle_as_generalised_element.h
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

oomph_add_executable(
  NAME old_macro_element_free_boundary_poisson_for_doc
  SOURCES old_for_doc.cc circle_as_generalised_element.h
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME interaction.free_boundary_poisson
  DEPENDS_ON algebraic_free_boundary_poisson
             circle
             doc_sparse_macro_node_update
             elastic_poisson
             free_boundary_poisson_with_elastic_mesh_update
             geom_object_element
             macro_element_free_boundary_poisson
             macro_element_free_boundary_poisson_non_ref
             old_macro_element_free_boundary_poisson_for_doc
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata)

message(VERBOSE "Leaving free_boundary_poisson subdirectory")
# ------------------------------------------------------------------------------
