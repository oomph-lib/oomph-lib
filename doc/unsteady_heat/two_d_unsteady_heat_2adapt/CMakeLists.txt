# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered two_d_unsteady_heat_2adapt subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(two_d_unsteady_heat_2adapt C CXX Fortran)
set(OOMPH_ROOT_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../")

include("${OOMPH_ROOT_DIR}/cmake/OomphGenerateDocFrom.cmake")
oomph_generate_doc_from(OOMPH_ROOT_DIR "${OOMPH_ROOT_DIR}"
                        DOCFILE two_d_unsteady_heat_2adapt.txt)

# The symlink we want to create
set(SRC_GIF
    "${CMAKE_CURRENT_SOURCE_DIR}/../two_d_unsteady_heat_ALE/figures/domain.gif")
set(DST_GIF "${CMAKE_CURRENT_SOURCE_DIR}/figures/domain.gif")

# A dummy stamp file so the rule has an OUTPUT CMake can track
set(STAMP "${CMAKE_CURRENT_BINARY_DIR}/domain_gif.stamp")

# Make sure we have a symlink
add_custom_command(
  OUTPUT "${STAMP}"
  COMMAND
    ${CMAKE_COMMAND} -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/ensure_symlink.cmake"
    "${SRC_GIF}" "${DST_GIF}"
  COMMAND ${CMAKE_COMMAND} -E touch "${STAMP}"
  DEPENDS "${SRC_GIF}"
  COMMENT "Ensuring symlink ${DST_GIF} --> ${SRC_GIF}"
  VERBATIM)

# Create a target that needs STAMP, which will trigger the above command
add_custom_target(generate_domain_gif ALL DEPENDS "${STAMP}")

message(VERBOSE "Leaving two_d_unsteady_heat_2adapt subdirectory")
# ------------------------------------------------------------------------------
