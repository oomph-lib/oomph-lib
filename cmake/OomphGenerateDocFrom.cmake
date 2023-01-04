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
include_guard()

# FIXME: Switch to using add_custom_command(...) instead of target to avoid
# constantly rebuilding...

# ------------------------------------------------------------------------------
function(oomph_generate_doc_from)
  # Define the supported set of keywords
  set(PREFIX ARG)
  set(FLAGS BUILD_DOCS_TARGET SUPPRESS_LATEX_IN_THIS_DIRECTORY)
  set(SINGLE_VALUE_ARGS OOMPH_ROOT_DIR DOCFILE)
  set(MULTI_VALUE_ARGS)

  # Process the arguments passed in
  include(CMakeParseArguments)
  cmake_parse_arguments(PARSE_ARGV 0 ${PREFIX} "${FLAGS}"
                        "${SINGLE_VALUE_ARGS}" "${MULTI_VALUE_ARGS}")

  # Strip the prefix
  set(OOMPH_ROOT_DIR ${${PREFIX}_OOMPH_ROOT_DIR})
  set(DOCFILE ${${PREFIX}_DOCFILE})
  set(BUILD_DOCS_TARGET ${${PREFIX}_BUILD_DOCS_TARGET})
  set(SUPPRESS_LATEX_IN_THIS_DIRECTORY
      ${${PREFIX}_SUPPRESS_LATEX_IN_THIS_DIRECTORY})

  # TODO: Add sanity checks to make sure the oomph-lib root is correct

  # Look for Doxygen (mandatory) and pdflatex (not mandatory)
  if(NOT Doxygen_FOUND)
    find_package(Doxygen 1.9.6 REQUIRED)
  endif()
  if(NOT DEFINED PATH_TO_PDFLATEX)
    find_program(PATH_TO_PDFLATEX NAMES pdflatex)
  endif()

  # NOTE: We can't verify that the docfile exists straight away as it might only
  # get generated at build-time; see e.g. doc/index/oomph_index.txt, but we can
  # check for it at build-time when it is needed.

  # Get relative path to oomph-lib root
  cmake_path(
    RELATIVE_PATH OOMPH_ROOT_DIR BASE_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    OUTPUT_VARIABLE RELATIVE_PATH_TO_OOMPH_ROOT_DIR)

  # Strip trailing slash (not necessary, just for prettiness...)
  string(REGEX REPLACE "/$" "" RELATIVE_PATH_TO_OOMPH_ROOT_DIR
                       "${RELATIVE_PATH_TO_OOMPH_ROOT_DIR}")

  # Extract the filename stem and extension; make sure the docfile is a text
  # file
  cmake_path(GET DOCFILE EXTENSION DOCFILE_EXT)
  cmake_path(GET DOCFILE STEM DOCFILE_STEM)
  if(NOT (DOCFILE_EXT MATCHES ".(txt|TXT)"))
    message(
      FATAL_ERROR
        "Expected argument '${DOCFILE}' to oomph_generate_doc_from(...) to be a .txt file, not a ${DOCFILE_EXT} file!"
    )
  endif()

  # Generate doxygen-ified header from ${DOCFILE}
  #
  # FIXME: txt2h_new.sh has a hard reliance on the system having 'awk'; need to
  # check for gawk, nawk or mawk if its missing

  # Add a link in the doxygen-ified header to the to-be-generated PDF doc
  if(NOT SUPPRESS_LATEX_IN_THIS_DIRECTORY)
    add_custom_command(
      OUTPUT "${CMAKE_CURRENT_SOURCE_DIR}/${DOCFILE}_doxygenified.h"
      COMMAND "${OOMPH_ROOT_DIR}/scripts/txt2h_new.sh" ${DOCFILE}
      COMMAND cp ${DOCFILE}_doxygenified.h ${DOCFILE}_doxygenified.h.junk
      COMMAND
        sed
        [=[s/\*\*\//\n<hr>\n<hr>\n\\\section pdf PDF file\nA <a href=\"..\/latex\/refman.pdf\">pdf version<\/a> of this document is available.\n\*\*\//g]=]
        ${DOCFILE}_doxygenified.h.junk > ${DOCFILE}_doxygenified.h
      COMMAND rm -f ${DOCFILE}_doxygenified.h.junk
      DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/${DOCFILE}"
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      VERBATIM)
  else()
    add_custom_command(
      OUTPUT "${CMAKE_CURRENT_SOURCE_DIR}/${DOCFILE}_doxygenified.h"
      COMMAND "${OOMPH_ROOT_DIR}/scripts/txt2h_new.sh" ${DOCFILE}
      DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/${DOCFILE}"
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      VERBATIM)
  endif()

  # Generate header and footer
  add_custom_command(
    OUTPUT "${CMAKE_CURRENT_SOURCE_DIR}/oomph-lib_header.html"
           "${CMAKE_CURRENT_SOURCE_DIR}/oomph-lib_footer.html"
    COMMAND "${OOMPH_ROOT_DIR}/scripts/build_oomph_html_header.sh"
            "${RELATIVE_PATH_TO_OOMPH_ROOT_DIR}"
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    VERBATIM)

  # Pull "doc/extra_latex_style_files/" folder into the current folder. If we're
  # in the doc/ folder, don't do anything
  set(LATEX_STYLE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/extra_latex_style_files")
  if((NOT EXISTS "${LATEX_STYLE_FILES}") OR (IS_SYMLINK "${LATEX_STYLE_FILES}"))
    add_custom_command(
      OUTPUT "${CMAKE_CURRENT_SOURCE_DIR}/extra_latex_style_files"
      COMMAND ln -sf "${OOMPH_ROOT_DIR}/doc/extra_latex_style_files"
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      VERBATIM)
  endif()

  # Hash the path to create a unique ID for our targets but shorten it to the
  # first 7 characters for brevity. A unique ID is required to avoid clashes
  # with targets in other directories
  string(SHA1 PATH_HASH "${CMAKE_CURRENT_SOURCE_DIR}")
  string(SUBSTRING ${PATH_HASH} 0 7 PATH_HASH)

  # After the doxygen-ified header and header/footer have been generated, start
  # the build of html/index.html
  add_custom_target(
    build_docs_${PATH_HASH} ALL
    DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/${DOCFILE}_doxygenified.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/oomph-lib_header.html"
            "${CMAKE_CURRENT_SOURCE_DIR}/oomph-lib_footer.html"
            "${CMAKE_CURRENT_SOURCE_DIR}/extra_latex_style_files")

  # Set the flag BUILD_DOCS_TARGET to the name of the build_docs_* target incase
  # the user needs to add their own dependencies
  if(BUILD_DOCS_TARGET)
    set(BUILD_DOCS_TARGET build_docs_${PATH_HASH} PARENT_SCOPE)
  endif()

  # Build! Build! Build!
  add_custom_command(
    TARGET build_docs_${PATH_HASH}
    POST_BUILD
    COMMAND doxygen
    COMMAND cp "${OOMPH_ROOT_DIR}/doc/figures/doxygen.png" html
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    BYPRODUCTS "${CMAKE_CURRENT_SOURCE_DIR}/html/"
               "${CMAKE_CURRENT_SOURCE_DIR}/latex/"
    VERBATIM)

  # To ensure we can clean up the html/ and latex/ directories (i.e. not plain
  # files), we have to set the ADDITIONAL_CLEAN_FILES property. Bit odd...
  set_property(
    TARGET build_docs_${PATH_HASH}
    APPEND
    PROPERTY ADDITIONAL_CLEAN_FILES "${CMAKE_CURRENT_SOURCE_DIR}/html/"
             "${CMAKE_CURRENT_SOURCE_DIR}/latex/")

  # If desired, build the PDF documentation
  if(PATH_TO_PDFLATEX AND (NOT SUPPRESS_LATEX_IN_THIS_DIRECTORY))
    add_custom_command(
      TARGET build_docs_${PATH_HASH}
      POST_BUILD
      COMMAND cd latex
      COMMAND rm -f refman.pdf
      COMMAND echo "\\\\end{document}" >> index.tex
      COMMAND mv refman.tex refman.tex.back
      COMMAND "${OOMPH_ROOT_DIR}/scripts/customise_latex.bash" refman.tex.back >
              refman.tex
      COMMAND "${OOMPH_ROOT_DIR}/scripts/tweak_doxygen_latex_style_file.bash"
      COMMAND make -i > /dev/null || make -i || echo "PDF build failed!"
      COMMAND mv refman.pdf ../${DOCFILE_STEM}.pdf
      COMMAND ln -s ../${DOCFILE_STEM}.pdf refman.pdf
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      BYPRODUCTS "${CMAKE_CURRENT_SOURCE_DIR}/${DOCFILE_STEM}.pdf"
      VERBATIM)
  endif()
endfunction()
# ------------------------------------------------------------------------------
