# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered two_d_unsteady_heat_2adapt subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(two_d_unsteady_heat_2adapt C CXX Fortran)
set(OOMPH_ROOT_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../")

include("${OOMPH_ROOT_DIR}/cmake/OomphGenerateDocFrom.cmake")
oomph_generate_doc_from(
  OOMPH_ROOT_DIR "${OOMPH_ROOT_DIR}"
  DOCFILE two_d_unsteady_heat_2adapt.txt
  DEFINE_BUILD_DOCS_TARGET_IN_CURRENT_SCOPE)

# The symlink we want to create
set(SRC_GIF
    "${CMAKE_CURRENT_SOURCE_DIR}/../two_d_unsteady_heat_ALE/figures/domain.gif")
set(DST_GIF "${CMAKE_CURRENT_SOURCE_DIR}/figures/domain.gif")

# Make sure we have a symlink
add_custom_command(
  OUTPUT "${DST_GIF}"
  COMMAND ${CMAKE_COMMAND} -E copy_if_different "${SRC_GIF}" "${DST_GIF}"
  DEPENDS "${SRC_GIF}"
  COMMENT "Copying domain.gif for documentation"
  VERBATIM)

# Create a custom target to drive the copy. When the target is built, it relies
# on the copied image existing, which will in turn trigger the above command
add_custom_target(copy_domain_gif ALL DEPENDS "${DST_GIF}")

# Make sure copy_domain_gif runs before we build the documentation
add_dependencies(${BUILD_DOCS_TARGET} copy_domain_gif)

message(VERBOSE "Leaving two_d_unsteady_heat_2adapt subdirectory")
# ------------------------------------------------------------------------------
