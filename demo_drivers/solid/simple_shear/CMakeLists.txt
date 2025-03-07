# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered simple_shear subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(simple_shear C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

# Precompile the mesh as an OBJECT library
add_library(precompiled_mesh OBJECT precompiled_mesh.cc)
target_link_libraries(precompiled_mesh PRIVATE oomph::solid)
target_compile_features(precompiled_mesh
                        PRIVATE cxx_std_${OOMPH_CMAKE_CXX_STANDARD})

oomph_add_executable(
  NAME simple_shear
  SOURCES simple_shear.cc
  LIBRARIES precompiled_mesh oomph::meshes oomph::generic oomph::solid
            oomph::constitutive)

oomph_add_executable(
  NAME refineable_simple_shear
  SOURCES refineable_simple_shear.cc
  LIBRARIES precompiled_mesh oomph::meshes oomph::generic oomph::solid
            oomph::constitutive)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME solid.simple_shear
  DEPENDS_ON simple_shear refineable_simple_shear
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata
  LABELS solid simple_shear)

# This test needs more than the default allotted time to complete
set_tests_properties(solid.simple_shear PROPERTIES TIMEOUT 3600)

message(VERBOSE "Leaving simple_shear subdirectory")
# ------------------------------------------------------------------------------
