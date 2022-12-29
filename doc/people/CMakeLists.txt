# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered people subdirectory")

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
project(people C CXX Fortran)
if(PROJECT_IS_TOP_LEVEL)
  find_package(oomphlib CONFIG REQUIRED PATHS "../../install")
endif()
find_package(Doxygen 1.9.2 REQUIRED)

set(DOCFILE people)

set(SUPPRESS_LATEX_IN_THIS_DIRECTORY FALSE)

find_program(PATH_TO_PDFLATEX NAMES pdflatex)

add_custom_command(
  OUTPUT ${DOCFILE}.txt_doxygenified.h
  COMMAND "${OOMPH_ROOT_DIR}/scripts/txt2h_new.sh" ${DOCFILE}.txt
  WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
  VERBATIM)

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
add_custom_target(
  build_docs ALL
  DEPENDS ${DOCFILE}.txt_doxygenified.h
          "${CMAKE_CURRENT_SOURCE_DIR}/oomph-lib_header.html"
          "${CMAKE_CURRENT_SOURCE_DIR}/oomph-lib_footer.html"
  COMMAND ln -sf "${OOMPH_ROOT_DIR}/doc/extra_latex_style_files"
  WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
  VERBATIM)

if(NOT SUPPRESS_LATEX_IN_THIS_DIRECTORY)
  add_custom_command(
    TARGET build_docs
    POST_BUILD
    COMMAND cp ${DOCFILE}.txt_doxygenified.h ${DOCFILE}.txt_doxygenified.h.junk
    COMMAND
      sed
      [=[s/\*\*\//\n<hr>\n<hr>\n\\\section pdf PDF file\nA <a href=\"..\/latex\/refman.pdf\">pdf version<\/a> of this document is available.\n\*\*\//g]=]
      ${DOCFILE}.txt_doxygenified.h.junk > ${DOCFILE}.txt_doxygenified.h
    COMMAND rm ${DOCFILE}.txt_doxygenified.h.junk
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    VERBATIM)
endif()

add_custom_command(
  TARGET build_docs
  POST_BUILD
  COMMAND doxygen
  COMMAND cp "${OOMPH_ROOT_DIR}/doc/figures/doxygen.png" html
  WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
  VERBATIM)

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
    COMMAND mv refman.pdf ../${DOCFILE}.pdf
    COMMAND ln -s ../${DOCFILE}.pdf refman.pdf
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    VERBATIM)
endif()

message(VERBOSE "Leaving people subdirectory")
# ------------------------------------------------------------------------------
