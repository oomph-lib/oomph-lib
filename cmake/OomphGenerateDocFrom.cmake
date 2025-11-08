include_guard()

function(oomph_generate_doc_from)
  # ----------------------- argument parsing -----------------------
  set(PREFIX ARG)
  set(FLAGS DEFINE_BUILD_DOCS_TARGET_IN_CURRENT_SCOPE
      SUPPRESS_LATEX_IN_THIS_DIRECTORY)
  set(SINGLE_VALUE_ARGS OOMPH_ROOT_DIR DOCFILE)
  include(CMakeParseArguments)
  cmake_parse_arguments(PARSE_ARGV 0 ${PREFIX} "${FLAGS}"
                        "${SINGLE_VALUE_ARGS}" "")

  set(OOMPH_ROOT_DIR "${ARG_OOMPH_ROOT_DIR}")
  set(DOCFILE "${ARG_DOCFILE}")
  set(SUPPRESS_PDF "${ARG_SUPPRESS_LATEX_IN_THIS_DIRECTORY}")
  set(DEFINE_BUILD_DOCS_TARGET_IN_CURRENT_SCOPE
      "${ARG_DEFINE_BUILD_DOCS_TARGET_IN_CURRENT_SCOPE}")

  find_package(Doxygen 1.9.6 REQUIRED)
  if(NOT DEFINED PATH_TO_PDFLATEX)
    find_program(PATH_TO_PDFLATEX NAMES pdflatex)
  endif()

  cmake_path(GET DOCFILE EXTENSION _ext)
  cmake_path(GET DOCFILE STEM DOC_STEM)
  if(NOT _ext MATCHES "^[.]?[tT][xX][tT]$")
    message(FATAL_ERROR "DOCFILE must be .txt (got ${DOCFILE}).")
  endif()

  # Files
  set(DOXIFY_OUT "${CMAKE_CURRENT_SOURCE_DIR}/${DOCFILE}_doxygenified.h")
  set(HTML_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/oomph-lib_header.html")
  set(HTML_FOOTER "${CMAKE_CURRENT_SOURCE_DIR}/oomph-lib_footer.html")
  set(HTML_INDEX "${CMAKE_CURRENT_SOURCE_DIR}/html/index.html")
  set(PDF_OUT "${CMAKE_CURRENT_SOURCE_DIR}/${DOC_STEM}.pdf")
  set(STYLE_LINK "${CMAKE_CURRENT_SOURCE_DIR}/extra_latex_style_files")

  cmake_path(RELATIVE_PATH OOMPH_ROOT_DIR BASE_DIRECTORY
             "${CMAKE_CURRENT_SOURCE_DIR}" OUTPUT_VARIABLE REL_PATH_TO_ROOT)
  string(REGEX REPLACE "/$" "" REL_PATH_TO_ROOT "${REL_PATH_TO_ROOT}")

  string(SHA1 PATH_HASH "${CMAKE_CURRENT_SOURCE_DIR}")
  string(SUBSTRING ${PATH_HASH} 0 7 PATH_HASH)

  if(SUPPRESS_PDF)
    add_custom_command(
      OUTPUT "${DOXIFY_OUT}"
      COMMAND "${OOMPH_ROOT_DIR}/scripts/txt2h_new.sh" "${DOCFILE}"
      DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/${DOCFILE}"
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      VERBATIM)
  else()
    add_custom_command(
      OUTPUT "${DOXIFY_OUT}"
      COMMAND "${OOMPH_ROOT_DIR}/scripts/txt2h_new.sh" "${DOCFILE}"
      COMMAND cp "${DOCFILE}_doxygenified.h" "${DOCFILE}_doxygenified.h.junk"
      COMMAND
        /bin/sh -c
        "sed -e 's@\\*\\*\\/@\\n<hr>\\n<hr>\\n\\\\section pdf PDF file\\nA <a href=\"../latex/refman.pdf\">pdf version</a> of this document is available.\\n\\\\*/@' \
               ${DOCFILE}_doxygenified.h.junk > ${DOCFILE}_doxygenified.h"
      COMMAND rm -f "${DOCFILE}_doxygenified.h.junk"
      DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/${DOCFILE}"
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      VERBATIM)
  endif()

  # HTML header / footer
  add_custom_command(
    OUTPUT "${HTML_HEADER}" "${HTML_FOOTER}"
    COMMAND "${OOMPH_ROOT_DIR}/scripts/build_oomph_html_header.sh"
            "${REL_PATH_TO_ROOT}"
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    VERBATIM)

  # The targets we'll build to generate the docs
  set(TARGETS_REQUIRED_FOR_BUILD_DOCS "${DOXIFY_OUT}" "${HTML_HEADER}"
      "${HTML_FOOTER}")

  # LaTeX style symlink
  if(NOT EXISTS "${STYLE_LINK}")
    add_custom_command(
      OUTPUT "${STYLE_LINK}"
      COMMAND ${CMAKE_COMMAND} -E create_symlink
              "${OOMPH_ROOT_DIR}/doc/extra_latex_style_files" "${STYLE_LINK}"
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      VERBATIM)

    # We're also generating a symlink ourselves
    list(APPEND TARGETS_REQUIRED_FOR_BUILD_DOCS "${STYLE_LINK}")
  endif()

  # Run Doxygen (when any of its deps have changed)
  add_custom_command(
    OUTPUT "${HTML_INDEX}"
    COMMAND doxygen -q
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${OOMPH_ROOT_DIR}/doc/figures/doxygen.png" html/doxygen.png
    DEPENDS ${TARGETS_REQUIRED_FOR_BUILD_DOCS}
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    BYPRODUCTS "${CMAKE_CURRENT_SOURCE_DIR}/html"
               "${CMAKE_CURRENT_SOURCE_DIR}/latex"
    COMMENT "Generating HTML & LaTeX documentation"
    VERBATIM)

  # We're always going to output the index.html
  set(DOC_OUTPUTS "${HTML_INDEX}")

  # Optional PDF build (depend on "html/index.html" to determine when to rerun)
  if(PATH_TO_PDFLATEX AND NOT SUPPRESS_PDF)
    add_custom_command(
      OUTPUT "${PDF_OUT}"
      COMMAND
        ${CMAKE_COMMAND} -E chdir latex /bin/sh -c
        " \
        rm -f refman.pdf && \
        echo '\\\\end{document}' >> index.tex && \
        mv refman.tex refman.tex.back && \
        '${OOMPH_ROOT_DIR}/scripts/customise_latex.bash' refman.tex.back > refman.tex && \
        '${OOMPH_ROOT_DIR}/scripts/tweak_doxygen_latex_style_file.bash' && \
        make -s -i \
             PDFLATEX='pdflatex -interaction=batchmode -halt-on-error' \
             MAKEINDEX='makeindex -q' \
          || { echo 'PDF build failed!  See latex/refman.log'; false; } && \
        cp refman.pdf '${PDF_OUT}' && \
        ln -sf '../${DOC_STEM}.pdf' refman.pdf \
      "
      DEPENDS "${HTML_INDEX}"
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      COMMENT "Building PDF documentation"
      VERBATIM)

    # ...and we'll also output a PDF
    list(APPEND DOC_OUTPUTS "${PDF_OUT}")
  endif()

  # Define the target that needs to be built
  add_custom_target(build_docs_${PATH_HASH} ALL DEPENDS ${DOC_OUTPUTS})

  # The "stuff" we have to make sure to clean up. Normally you can tell CMake
  # what the byproducts are and it will help clean it up, but these are
  # directories and CMake objects to deleting non-empty directories with the
  # typical workflow. To get around this, we will tell it to make sure to clean
  # up these extra folders (and possibly symlinks)
  set(ADDITIONAL_CLEAN_UP "${CMAKE_CURRENT_SOURCE_DIR}/html"
      "${CMAKE_CURRENT_SOURCE_DIR}/latex")

  # If we're also generating a symlink
  if(NOT EXISTS "${STYLE_LINK}")
    list(APPEND ADDITIONAL_CLEAN_UP "${STYLE_LINK}")
  endif()

  # Now tell CMake to make sure to clean these files up
  set_property(
    TARGET build_docs_${PATH_HASH}
    APPEND
    PROPERTY ADDITIONAL_CLEAN_FILES ${ADDITIONAL_CLEAN_UP})

  # Make the target name visible to the outside scope
  if(DEFINE_BUILD_DOCS_TARGET_IN_CURRENT_SCOPE)
    set(BUILD_DOCS_TARGET "build_docs_${PATH_HASH}" PARENT_SCOPE)
  endif()
endfunction()
