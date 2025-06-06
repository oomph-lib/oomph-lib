include_guard()

function(oomph_generate_doc_from)
  # ----------------------- argument parsing -----------------------
  set(PREFIX ARG)
  set(FLAGS BUILD_DOCS_TARGET SUPPRESS_LATEX_IN_THIS_DIRECTORY)
  set(SINGLE_VALUE_ARGS OOMPH_ROOT_DIR DOCFILE)
  include(CMakeParseArguments)
  cmake_parse_arguments(PARSE_ARGV 0 ${PREFIX} "${FLAGS}"
                        "${SINGLE_VALUE_ARGS}" "")

  set(OOMPH_ROOT_DIR "${ARG_OOMPH_ROOT_DIR}")
  set(DOCFILE "${ARG_DOCFILE}")
  set(SUPPRESS_PDF "${ARG_SUPPRESS_LATEX_IN_THIS_DIRECTORY}")

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

  # LaTeX style symlink
  if(NOT EXISTS "${STYLE_LINK}" AND NOT IS_SYMLINK "${STYLE_LINK}")
    add_custom_command(
      OUTPUT "${STYLE_LINK}"
      COMMAND ${CMAKE_COMMAND} -E create_symlink
              "${OOMPH_ROOT_DIR}/doc/extra_latex_style_files" "${STYLE_LINK}"
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      VERBATIM)
  endif()

  # Run doxygen
  add_custom_target(
    build_docs_${PATH_HASH} ALL DEPENDS "${DOXIFY_OUT}" "${HTML_HEADER}"
                                        "${HTML_FOOTER}" "${STYLE_LINK}")

  add_custom_command(
    TARGET build_docs_${PATH_HASH}
    POST_BUILD
    COMMAND doxygen -q
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${OOMPH_ROOT_DIR}/doc/figures/doxygen.png" html/doxygen.png
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    BYPRODUCTS "${CMAKE_CURRENT_SOURCE_DIR}/html"
               "${CMAKE_CURRENT_SOURCE_DIR}/latex" OUTPUT_QUIET
    VERBATIM)

  # Optional PDF build
  if(PATH_TO_PDFLATEX AND NOT SUPPRESS_PDF)
    add_custom_command(
      TARGET build_docs_${PATH_HASH}
      POST_BUILD
      COMMAND
        ${CMAKE_COMMAND} -E chdir latex /bin/sh -c
        "\
        rm -f refman.pdf &&\
        echo '\\\\end{document}' >> index.tex &&\
        mv refman.tex refman.tex.back &&\
        '${OOMPH_ROOT_DIR}/scripts/customise_latex.bash' refman.tex.back > refman.tex &&\
        '${OOMPH_ROOT_DIR}/scripts/tweak_doxygen_latex_style_file.bash' && \
        make -s -i \
             PDFLATEX='pdflatex -interaction=batchmode -halt-on-error' \
             MAKEINDEX='makeindex -q' \
          || { echo 'PDF build failed!  See latex/refman.log'; false; } &&\
        cp refman.pdf '${PDF_OUT}' &&\
        ln -sf '../${DOC_STEM}.pdf' refman.pdf\
      "
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      BYPRODUCTS "${PDF_OUT}"
      VERBATIM)
  endif()

  # ------------------ F) allow clean to wipe dirs ----------------
  set_property(
    TARGET build_docs_${PATH_HASH}
    APPEND
    PROPERTY ADDITIONAL_CLEAN_FILES "${CMAKE_CURRENT_SOURCE_DIR}/html"
             "${CMAKE_CURRENT_SOURCE_DIR}/latex")

  # ------------------ G) return target name (flag API) ------------
  if(ARG_BUILD_DOCS_TARGET)
    set(BUILD_DOCS_TARGET "build_docs_${PATH_HASH}" PARENT_SCOPE)
  endif()
endfunction()
