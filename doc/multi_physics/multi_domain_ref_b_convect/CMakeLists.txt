# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered multi_domain_ref_b_convect subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(multi_domain_ref_b_convect C CXX Fortran)
set(OOMPH_ROOT_DIR "${CMAKE_CURRENT_LIST_DIR}/../../..")

include("${OOMPH_ROOT_DIR}/cmake/OomphGenerateDocFrom.cmake")
oomph_generate_doc_from(OOMPH_ROOT_DIR "${OOMPH_ROOT_DIR}"
                        DOCFILE multi_domain_ref_b_convect.txt)

message(VERBOSE "Leaving multi_domain_ref_b_convect subdirectory")
# ------------------------------------------------------------------------------
