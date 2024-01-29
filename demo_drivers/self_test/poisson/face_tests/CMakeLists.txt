# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered face_tests subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(face_tests C CXX Fortran)
if(NOT oomphlib_FOUND)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../../install")
endif()
include(CTest)

oomph_add_executable(
  NAME q_faces_3d
  SOURCES q_faces_3d.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

oomph_add_executable(
  NAME q_faces_2d
  SOURCES q_faces_2d.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

oomph_add_executable(
  NAME t_faces_2d
  SOURCES t_faces_2d.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

oomph_add_executable(
  NAME t_faces_3d
  SOURCES t_faces_3d.cc
  LIBRARIES oomph::poisson oomph::meshes oomph::generic)

# ---------------------------------- TESTING ----------------------------------

oomph_add_test(
  TEST_NAME self_test.poisson.face_tests
  DEPENDS_ON q_faces_3d q_faces_2d t_faces_2d t_faces_3d
  COMMAND ./validate.sh ${OOMPH_ROOT_DIR}
  TEST_FILES validate.sh validata
  LABELS face_tests poisson)

message(VERBOSE "Leaving face_tests subdirectory")
# ------------------------------------------------------------------------------
