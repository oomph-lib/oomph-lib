# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered eighth_sphere_poisson subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(eighth_sphere_poisson C CXX Fortran)
set(OOMPH_ROOT_DIR "${CMAKE_CURRENT_LIST_DIR}/../../..")

include("${OOMPH_ROOT_DIR}/cmake/OomphGenerateDocFrom.cmake")
oomph_generate_doc_from(OOMPH_ROOT_DIR "${OOMPH_ROOT_DIR}"
                        DOCFILE eighth_sphere_poisson.txt)

message(VERBOSE "Leaving eighth_sphere_poisson subdirectory")
# ------------------------------------------------------------------------------
