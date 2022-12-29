# cmake-format: off
# =============================================================================
# DESCRIPTION:
# ------------
#
#
# USAGE:
# ------
#     include(OomphGenerateDocFrom)
#     oomph_generate_doc_from(<docfile>)
#
# =============================================================================
# cmake-format: on

# ------------------------------------------------------------------------------
function(oomph_generate_doc_from)
  # Define the supported set of keywords
  set(PREFIX ARG)
  set(FLAGS)
  set(SINGLE_VALUE_ARGS OOMPH_ROOT_DIR DOCFILE)
  set(MULTI_VALUE_ARGS)

  # Process the arguments passed in
  include(CMakeParseArguments)
  cmake_parse_arguments(PARSE_ARGV 0 ${PREFIX} "${FLAGS}"
                        "${SINGLE_VALUE_ARGS}" "${MULTI_VALUE_ARGS}")

  # Strip the prefix
  set(OOMPH_ROOT_DIR ${${PREFIX}_OOMPH_ROOT_DIR})
  set(DOCFILE ${${PREFIX}_DOCFILE})

  # TODO: Add sanity checks to make sure the oomph-lib root is correct

  # Look for Doxygen
  if(NOT Doxygen_FOUND)
    find_package(Doxygen 1.9.2 REQUIRED)
  endif()

  # Extract the filename stem and extension; make sure it's a text file
  cmake_path(GET DOCFILE EXTENSION DOCFILE_EXT)
  cmake_path(GET DOCFILE STEM DOCFILE_STEM)
  if(NOT (DOCFILE_EXT MATCHES ".(txt|TXT)"))
    message(
      FATAL_ERROR
        "Expected argument '${DOCFILE}' to oomph_generate_doc_from(...) to be a .txt file, not a ${DOCFILE_EXT} file!"
    )
  endif()

  if(NOT DEFINED PATH_TO_PDFLATEX)
    find_program(PATH_TO_PDFLATEX NAMES pdflatex)
  endif()

  # Generate doxygen-ified header from ${DOCFILE_STEM}.txt
  add_custom_command(
    OUTPUT "${CMAKE_CURRENT_SOURCE_DIR}/${DOCFILE_STEM}.txt_doxygenified.h"
    COMMAND "${OOMPH_ROOT_DIR}/scripts/txt2h_new.sh" ${DOCFILE_STEM}.txt
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    VERBATIM)

  # Generate header and footer
  add_custom_command(
    OUTPUT "${CMAKE_CURRENT_SOURCE_DIR}/oomph-lib_header.html"
           "${CMAKE_CURRENT_SOURCE_DIR}/oomph-lib_footer.html"
    COMMAND "${OOMPH_ROOT_DIR}/scripts/build_oomph_html_header.sh"
            "${OOMPH_ROOT_DIR}"
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    VERBATIM)

  # FIXME: Handle the case that build_docs is already defined
  #
  # if(NOT TARGET build_docs) add_custom_target(build_docs ALL) endif()

  # After the doxygen-ified header and header/footer have been generated, start
  # the build of html/index.html
  add_custom_target(
    build_docs ALL
    DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/${DOCFILE_STEM}.txt_doxygenified.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/oomph-lib_header.html"
            "${CMAKE_CURRENT_SOURCE_DIR}/oomph-lib_footer.html"
    COMMAND ln -sf "${OOMPH_ROOT_DIR}/doc/extra_latex_style_files"
    BYPRODUCTS "${CMAKE_CURRENT_SOURCE_DIR}/extra_latex_style_files"
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    VERBATIM)

  # Add a link in the doxygen-ified header to the to-be-generated PDF doc
  if(NOT SUPPRESS_LATEX_IN_THIS_DIRECTORY)
    add_custom_command(
      TARGET build_docs
      POST_BUILD
      COMMAND cp ${DOCFILE_STEM}.txt_doxygenified.h
              ${DOCFILE_STEM}.txt_doxygenified.h.junk
      COMMAND
        sed
        [=[s/\*\*\//\n<hr>\n<hr>\n\\\section pdf PDF file\nA <a href=\"..\/latex\/refman.pdf\">pdf version<\/a> of this document is available.\n\*\*\//g]=]
        ${DOCFILE_STEM}.txt_doxygenified.h.junk >
        ${DOCFILE_STEM}.txt_doxygenified.h
      COMMAND rm ${DOCFILE_STEM}.txt_doxygenified.h.junk
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      VERBATIM)
  endif()

  # Build! Build! Build!
  add_custom_command(
    TARGET build_docs
    POST_BUILD
    COMMAND doxygen
    COMMAND cp "${OOMPH_ROOT_DIR}/doc/figures/doxygen.png" html
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    BYPRODUCTS "${CMAKE_CURRENT_SOURCE_DIR}/html/"
               "${CMAKE_CURRENT_SOURCE_DIR}/latex/"
    VERBATIM)

  set_property(
    TARGET build_docs
    APPEND
    PROPERTY ADDITIONAL_CLEAN_FILES "${CMAKE_CURRENT_SOURCE_DIR}/html/"
             "${CMAKE_CURRENT_SOURCE_DIR}/latex/")

  # If desired, build the PDF documentation
  if(PATH_TO_PDFLATEX AND (NOT SUPPRESS_LATEX_IN_THIS_DIRECTORY))
    add_custom_command(
      TARGET build_docs
      POST_BUILD
      COMMAND cd latex
      COMMAND rm -f refman.pdf
      COMMAND echo "\\end{document}" >> index.tex
      COMMAND mv refman.tex refman.tex.back
      COMMAND "${OOMPH_ROOT_DIR}/scripts/customise_latex.bash" refman.tex.back >
              refman.tex
      COMMAND "${OOMPH_ROOT_DIR}/scripts/tweak_doxygen_latex_style_file.bash"
      COMMAND make -i > /dev/null || make -i
      COMMAND mv refman.pdf ../${DOCFILE_STEM}.pdf
      COMMAND ln -s ../${DOCFILE_STEM}.pdf refman.pdf
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      BYPRODUCTS "${CMAKE_CURRENT_SOURCE_DIR}/${DOCFILE_STEM}.pdf"
      VERBATIM)
  endif()

endfunction()
# ------------------------------------------------------------------------------
