# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered two_d_unsteady_heat_2adapt subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(two_d_unsteady_heat_2adapt C CXX Fortran)
set(OOMPH_ROOT_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../")

include("${OOMPH_ROOT_DIR}/cmake/OomphGenerateDocFrom.cmake")
oomph_generate_doc_from(OOMPH_ROOT_DIR "${OOMPH_ROOT_DIR}"
                        DOCFILE two_d_unsteady_heat_2adapt.txt)

add_custom_target(
  generate_domain_gif ALL
  COMMAND ${CMAKE_COMMAND} -E remove_if_exists
          "${CMAKE_CURRENT_SOURCE_DIR}/figures/domain.gif"
  COMMAND
    ${CMAKE_COMMAND} -E create_symlink
    "${CMAKE_CURRENT_SOURCE_DIR}/../two_d_unsteady_heat_ALE/figures/domain.gif"
    "${CMAKE_CURRENT_SOURCE_DIR}/figures/domain.gif"
  WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")

message(VERBOSE "Leaving two_d_unsteady_heat_2adapt subdirectory")
# ------------------------------------------------------------------------------
