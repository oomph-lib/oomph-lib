# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered third_party_meshes subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(third_party_meshes C CXX Fortran)
set(OOMPH_ROOT_DIR "${CMAKE_CURRENT_LIST_DIR}/../../..")

include("${OOMPH_ROOT_DIR}/cmake/OomphGenerateDocFrom.cmake")
oomph_generate_doc_from(OOMPH_ROOT_DIR "${OOMPH_ROOT_DIR}"
                        DOCFILE third_party_meshes.txt)

message(VERBOSE "Leaving third_party_meshes subdirectory")
# ------------------------------------------------------------------------------
