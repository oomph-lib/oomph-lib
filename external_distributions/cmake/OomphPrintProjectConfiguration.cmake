# cmake-format: off
# =============================================================================
# DESCRIPTION:
# ------------
# Prints the project configuration variables and compiler definitions.
#
# USAGE:
# ------
#   include(OomphPrintProjectConfiguration)
#   oomph_print_project_configuration()
# =============================================================================
include_guard()

# ------------------------------------------------------------------------------
function(oomph_print_project_configuration)
  # Define the supported set of keywords
  set(PREFIX ARG)
  set(FLAGS ENABLE_ALSO_PRINT_SETTINGS_AFTER_INSTALL ENABLE_SAVE_TO_FILE)
  set(SINGLE_VALUE_ARGS)
  set(MULTI_VALUE_ARGS)
  cmake_parse_arguments(PARSE_ARGV 0 ${PREFIX} "${FLAGS}" "${SINGLE_VALUE_ARGS}" "${MULTI_VALUE_ARGS}")

  # Combine the literal and the list variable into one list.
  list(PREPEND OOMPH_CONFIG_VARS CMAKE_BUILD_TYPE)

  # Colourising
  if(NOT WIN32)
    string(ASCII 27 Esc)
    set(RESET "${Esc}[m")
    set(BOLD "${Esc}[1m")
    set(RED "${Esc}[31m")
    set(GREEN "${Esc}[32m")
    set(YELLOW "${Esc}[33m")
    set(BLUE "${Esc}[34m")
    set(MAGENTA "${Esc}[35m")
    set(CYAN "${Esc}[36m")
    set(WHITE "${Esc}[37m")
    set(BOLD_RED "${Esc}[1;31m")
    set(BOLD_GREEN "${Esc}[1;32m")
    set(BOLD_YELLOW "${Esc}[1;33m")
    set(BOLD_BLUE "${Esc}[1;34m")
    set(BOLD_MAGENTA "${Esc}[1;35m")
    set(BOLD_CYAN "${Esc}[1;36m")
    set(BOLD_WHITE "${Esc}[1;37m")
  endif()

  # Initialise
  set(MARKER "⦿ ")
  set(OOMPH_SETTINGS_MESSAGE "\n")

  # Append configuration options
  string(
    APPEND
    OOMPH_SETTINGS_MESSAGE
    "********************************************************************************\n"
  )
  string(APPEND OOMPH_SETTINGS_MESSAGE
         "${BOLD}OOMPH-LIB THIRD-PARTY LIBRARIES OPTIONS:${RESET}\n")
  string(
    APPEND
    OOMPH_SETTINGS_MESSAGE
    "********************************************************************************\n"
  )
  foreach(OPTION IN LISTS OOMPH_CONFIG_VARS)
    if(NOT ${OPTION})
      string(
        APPEND OOMPH_SETTINGS_MESSAGE
        "  ${RED}${MARKER}${RESET} ${OPTION}: ${RED}'${${OPTION}}'${RESET}\n")
    else()
      string(
        APPEND
        OOMPH_SETTINGS_MESSAGE
        "  ${GREEN}${MARKER}${RESET} ${OPTION}: ${GREEN}'${${OPTION}}'${RESET}\n"
      )
    endif()
  endforeach()
  string(
    APPEND
    OOMPH_SETTINGS_MESSAGE
    "********************************************************************************\n"
  )

  # Print it
  file(RELATIVE_PATH OOMPH_RELPATH_TO_TPL_BUILD_DIR ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
  message(NOTICE "${OOMPH_SETTINGS_MESSAGE}")
  message(
    STATUS
      "${BOLD_MAGENTA}Project configured! Don't forget to run the build step with:\n\n"
      "\tcmake --build ${OOMPH_RELPATH_TO_TPL_BUILD_DIR}${RESET}\n")

  # Log to file if needed
  set(OOMPH_TPL_CONFIG_FILE "${CMAKE_CURRENT_BINARY_DIR}/oomph-lib-third-party-libraries-config.log")
  if(${PREFIX}_ENABLE_SAVE_TO_FILE OR ${PREFIX}_ENABLE_ALSO_PRINT_SETTINGS_AFTER_INSTALL)
    file(WRITE "${OOMPH_TPL_CONFIG_FILE}" "${OOMPH_SETTINGS_MESSAGE}")
  endif()
endfunction()
# ------------------------------------------------------------------------------
# cmake-format: on
