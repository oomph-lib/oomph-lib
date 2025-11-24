# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered fish_poisson2 subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(fish_poisson2 C CXX Fortran)
set(OOMPH_ROOT_DIR "${CMAKE_CURRENT_LIST_DIR}/../../..")

include("${OOMPH_ROOT_DIR}/cmake/OomphGenerateDocFrom.cmake")
oomph_generate_doc_from(OOMPH_ROOT_DIR "${OOMPH_ROOT_DIR}"
                        DOCFILE fish_poisson2.txt)

message(VERBOSE "Leaving fish_poisson2 subdirectory")
# ------------------------------------------------------------------------------
