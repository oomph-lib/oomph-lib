# cmake-format: off
# =============================================================================
# DESCRIPTION:
# ------------
# Prints the project configuration variables and compiler definitions.
#
# USAGE:
# ------
#          include(OomphPrintProjectSettings)
#          oomph_print_project_settings()
# =============================================================================
# cmake-format: on

# ------------------------------------------------------------------------------
function(oomph_print_project_settings)
  # Define the supported set of keywords
  set(PREFIX ARG)
  set(FLAGS ENABLE_PRINT_AFTER_INSTALL ENABLE_SAVE_TO_FILE)
  set(SINGLE_VALUE_ARGS)
  set(MULTI_VALUE_ARGS)
  cmake_parse_arguments(PARSE_ARGV 0 ${PREFIX} "${FLAGS}"
                        "${SINGLE_VALUE_ARGS}" "${MULTI_VALUE_ARGS}")

  # Initialise
  set(MARKER "⦿")
  set(OOMPH_SETTINGS_MESSAGE "\n")

  # Append configuration options
  string(APPEND OOMPH_SETTINGS_MESSAGE
         "************************************************************\n")
  string(APPEND OOMPH_SETTINGS_MESSAGE "================================\n")
  string(APPEND OOMPH_SETTINGS_MESSAGE "OOMPH-LIB CONFIGURATION OPTIONS:\n")
  string(APPEND OOMPH_SETTINGS_MESSAGE "================================\n")
  foreach(OPTION IN LISTS OOMPH_CONFIG_VARS)
    string(APPEND OOMPH_SETTINGS_MESSAGE
           "  ${MARKER} ${OPTION}: '${${OPTION}}'\n")
  endforeach()

  # Get the list compile definitions
  get_directory_property(
    OOMPH_COMPILE_DEFINITIONS DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                                        COMPILE_DEFINITIONS)

  # Sort the list alphabetically
  list(SORT OOMPH_COMPILE_DEFINITIONS)

  # Now print the compile definitions
  string(APPEND OOMPH_SETTINGS_MESSAGE "\n")
  string(APPEND OOMPH_SETTINGS_MESSAGE "===============================\n")
  string(APPEND OOMPH_SETTINGS_MESSAGE "OOMPH-LIB COMPILER DEFINITIONS:\n")
  string(APPEND OOMPH_SETTINGS_MESSAGE "===============================\n")
  foreach(DEFN ${OOMPH_COMPILE_DEFINITIONS})
    string(APPEND OOMPH_SETTINGS_MESSAGE "  ${MARKER} ${DEFN}\n")
  endforeach()
  string(APPEND OOMPH_SETTINGS_MESSAGE
         "************************************************************\n")

  # Print it
  message(NOTICE "${OOMPH_SETTINGS_MESSAGE}")

  # Log to file if needed
  set(OUTPUT_FILE "${CMAKE_CURRENT_BINARY_DIR}/oomphlib-configuration.log")
  if(${PREFIX}_ENABLE_SAVE_TO_FILE OR ${PREFIX}_ENABLE_PRINT_AFTER_INSTALL)
    file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/oomphlib-configuration.log"
         "${OOMPH_SETTINGS_MESSAGE}")
  endif()

  # Print during install step; oomph_print_project_settings(...) must be called
  # at the end of the configuration step for this to get printed at the end
  if(${PREFIX}_ENABLE_PRINT_AFTER_INSTALL)
    install(CODE "execute_process(COMMAND cat \"${OUTPUT_FILE}\")")
  endif()
endfunction()
# ------------------------------------------------------------------------------
