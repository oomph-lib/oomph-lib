# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE
        "Entered distributed_linear_algebra_infrastructure subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(distributed_linear_algebra_infrastructure C CXX Fortran)
set(OOMPH_ROOT_DIR "${CMAKE_CURRENT_LIST_DIR}/../../..")

include("${OOMPH_ROOT_DIR}/cmake/OomphGenerateDocFrom.cmake")
oomph_generate_doc_from(OOMPH_ROOT_DIR "${OOMPH_ROOT_DIR}"
                        DOCFILE distributed_linear_algebra_infrastructure.txt)
message(VERBOSE
        "Leaving distributed_linear_algebra_infrastructure subdirectory")
# ------------------------------------------------------------------------------
