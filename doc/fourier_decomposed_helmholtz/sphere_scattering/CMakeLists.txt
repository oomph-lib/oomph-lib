# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered sphere_scattering subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(sphere_scattering C CXX Fortran)
set(OOMPH_ROOT_DIR "${CMAKE_CURRENT_LIST_DIR}/../../..")

include("${OOMPH_ROOT_DIR}/cmake/OomphGenerateDocFrom.cmake")
oomph_generate_doc_from(OOMPH_ROOT_DIR "${OOMPH_ROOT_DIR}"
                        DOCFILE sphere_scattering.txt)

message(VERBOSE "Leaving sphere_scattering subdirectory")
# ------------------------------------------------------------------------------
