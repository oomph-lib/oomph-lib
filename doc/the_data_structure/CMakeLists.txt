# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered the_data_structure subdirectory")

cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(the_data_structure C CXX Fortran)
set(OOMPH_ROOT_DIR "${CMAKE_CURRENT_LIST_DIR}/../../")

include("${OOMPH_ROOT_DIR}/cmake/OomphGenerateDocFrom.cmake")
oomph_generate_doc_from(
  OOMPH_ROOT_DIR "${OOMPH_ROOT_DIR}"
  DOCFILE the_data_structure.txt
  DEFINE_BUILD_DOCS_TARGET_IN_CURRENT_SCOPE)

# Patch html/ files after the docs have been built
add_custom_command(
  TARGET ${BUILD_DOCS_TARGET}
  POST_BUILD
  # Patch the class list (missing linebreak)
  COMMAND cp html/classes.html html/classes.html.back
  COMMAND
    sed
    [=[s/<div class="classindex">/<p style="color:White;">.<\/p><p style="color:White;">.<\/p><div class="classindex">/g]=]
    html/classes.html.back > html/classes.html
  COMMAND rm -f html/classes.html.back
  # Add the alphabetical index of member functions to the top of the list for
  # those starting with "a" and use that as the general landing page.
  COMMAND
    sed
    [=[s/they belong to:<\/div>/they belong to:<\/div><ul><li> <a href="functions_a.html">A<\/a><li> <a href="functions_b.html">B<\/a><li> <a href="functions_c.html">C<\/a><li> <a href="functions_d.html">D<\/a><li> <a href="functions_e.html">E<\/a><li> <a href="functions_f.html">F<\/a><li> <a href="functions_g.html">G<\/a><li> <a href="functions_h.html">H<\/a><li> <a href="functions_i.html">I<\/a><li> <a href="functions_j.html">J<\/a><li> <a href="functions_k.html">K<\/a><li> <a href="functions_l.html">L<\/a><li> <a href="functions_m.html">M<\/a><li> <a href="functions_n.html">N<\/a><li> <a href="functions_o.html">O<\/a><li> <a href="functions_p.html">P<\/a><li> <a href="functions_q.html">Q<\/a><li> <a href="functions_r.html">R<\/a><li> <a href="functions_s.html">S<\/a><li> <a href="functions_t.html">T<\/a><li> <a href="functions_u.html">U<\/a><li> <a href="functions_v.html">V<\/a><li> <a href="functions_x.html">X<\/a><li> <a href="functions_y.html">Y<\/a><li> <a href="functions_z.html">Z<\/a><\/ul>/g]=]
    html/functions_a.html > html/functions.html
  WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
  VERBATIM)

message(VERBOSE "Leaving the_data_structure subdirectory")
# ------------------------------------------------------------------------------
