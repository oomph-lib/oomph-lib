# cmake-format: off
# =============================================================================
# DESCRIPTION:
# ------------
# Docs the project configuration variables and compiler definitions.
#
# USAGE:
# ------
#   include(OomphDocProjectSettings)
#   oomph_doc_project_settings()
# =============================================================================
# cmake-format: on
include_guard()

# ------------------------------------------------------------------------------
function(oomph_doc_project_settings)
  # Define the supported set of keywords
  set(PREFIX ARG)
  set(FLAGS ENABLE_ALSO_PRINT_SETTINGS_AFTER_INSTALL ENABLE_SAVE_TO_FILE)
  set(SINGLE_VALUE_ARGS)
  set(MULTI_VALUE_ARGS)
  cmake_parse_arguments(PARSE_ARGV 0 ${PREFIX} "${FLAGS}"
                        "${SINGLE_VALUE_ARGS}" "${MULTI_VALUE_ARGS}")

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

  # Doc it
  message(NOTICE "${OOMPH_SETTINGS_MESSAGE}")

  # Log to file if needed
  set(OUTPUT_FILE
      "${CMAKE_CURRENT_BINARY_DIR}/oomphlib-third-party-libraries-info.log")
  if(${PREFIX}_ENABLE_SAVE_TO_FILE
     OR ${PREFIX}_ENABLE_ALSO_PRINT_SETTINGS_AFTER_INSTALL)
    file(WRITE
         "${CMAKE_CURRENT_BINARY_DIR}/oomphlib-third-party-libraries-info.log"
         "${OOMPH_SETTINGS_MESSAGE}")
  endif()

  # Doc during install step; oomph_doc_project_settings(...) must be called at
  # the end of the configuration step for this to get doced at the end
  if(${PREFIX}_ENABLE_ALSO_PRINT_SETTINGS_AFTER_INSTALL)
    install(CODE "execute_process(COMMAND cat \"${OUTPUT_FILE}\")")
  endif()
endfunction()
# ------------------------------------------------------------------------------
