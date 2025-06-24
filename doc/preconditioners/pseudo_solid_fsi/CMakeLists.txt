# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered pseudo_solid_fsi subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(pseudo_solid_fsi C CXX Fortran)
set(OOMPH_ROOT_DIR "${CMAKE_CURRENT_LIST_DIR}/../../..")

include("${OOMPH_ROOT_DIR}/cmake/OomphGenerateDocFrom.cmake")
oomph_generate_doc_from(OOMPH_ROOT_DIR "${OOMPH_ROOT_DIR}"
                        DOCFILE pseudo_solid_fsi.txt)

message(VERBOSE "Leaving pseudo_solid_fsi subdirectory")
# ------------------------------------------------------------------------------
