# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered pml_time_harmonic_linear_elasticity subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(pml_time_harmonic_linear_elasticity C CXX Fortran)
set(OOMPH_ROOT_DIR "${CMAKE_CURRENT_LIST_DIR}/../..")

include("${OOMPH_ROOT_DIR}/cmake/OomphGenerateDocFrom.cmake")
oomph_generate_doc_from(OOMPH_ROOT_DIR "${OOMPH_ROOT_DIR}"
                        DOCFILE pml_time_harmonic_linear_elasticity.txt)

message(VERBOSE "Leaving pml_time_harmonic_linear_elasticity subdirectory")
# ------------------------------------------------------------------------------
