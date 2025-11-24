# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered single_layer_free_surface subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(single_layer_free_surface C CXX Fortran)
set(OOMPH_ROOT_DIR "${CMAKE_CURRENT_LIST_DIR}/../../..")

include("${OOMPH_ROOT_DIR}/cmake/OomphGenerateDocFrom.cmake")
oomph_generate_doc_from(OOMPH_ROOT_DIR "${OOMPH_ROOT_DIR}"
                        DOCFILE single_layer_free_surface.txt)

message(VERBOSE "Leaving single_layer_free_surface subdirectory")
# ------------------------------------------------------------------------------
