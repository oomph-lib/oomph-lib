# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered index subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(index C CXX Fortran)

set(OOMPH_ROOT_DIR "${CMAKE_CURRENT_LIST_DIR}/../../")
include("${OOMPH_ROOT_DIR}/cmake/OomphGenerateDocFrom.cmake")
oomph_generate_doc_from(
  OOMPH_ROOT_DIR "${OOMPH_ROOT_DIR}"
  DOCFILE oomph_index.txt
  BUILD_DOXY_HEADER_TARGET)

# Make sure to construct oomph_index.txt from intro.txt before we try to create
# the doxygen-ified header from the docfile
add_custom_target(
  generate_oomph_index
  COMMAND cp intro.txt oomph_index.txt
  COMMAND "${OOMPH_ROOT_DIR}/scripts/make_index.py" oomph.idx
          "${OOMPH_ROOT_DIR}" >> oomph_index.txt
  BYPRODUCTS "${CMAKE_CURRENT_SOURCE_DIR}/oomph_index.txt"
  WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
add_dependencies(${BUILD_DOXY_HEADER_TARGET} generate_oomph_index)

message(VERBOSE "Leaving index subdirectory")
# ------------------------------------------------------------------------------
