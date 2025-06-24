# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered jeffery_orbit subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(jeffery_orbit C CXX Fortran)
set(OOMPH_ROOT_DIR "${CMAKE_CURRENT_LIST_DIR}/../../..")

include("${OOMPH_ROOT_DIR}/cmake/OomphGenerateDocFrom.cmake")
oomph_generate_doc_from(OOMPH_ROOT_DIR "${OOMPH_ROOT_DIR}"
                        DOCFILE jeffery_orbit.txt)

message(VERBOSE "Leaving jeffery_orbit subdirectory")
# ------------------------------------------------------------------------------
