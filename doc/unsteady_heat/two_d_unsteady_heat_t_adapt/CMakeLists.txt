# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered two_d_unsteady_heat_t_adapt subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(two_d_unsteady_heat_t_adapt C CXX Fortran)
set(OOMPH_ROOT_DIR "${CMAKE_CURRENT_LIST_DIR}/../../..")

include("${OOMPH_ROOT_DIR}/cmake/OomphGenerateDocFrom.cmake")
oomph_generate_doc_from(OOMPH_ROOT_DIR "${OOMPH_ROOT_DIR}"
                        DOCFILE two_d_unsteady_heat_t_adapt.txt)

message(VERBOSE "Leaving two_d_unsteady_heat_t_adapt subdirectory")
# ------------------------------------------------------------------------------
