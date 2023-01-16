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
  set(MARKER "⦿")
  include(FeatureSummary)
  message(NOTICE "")
  message(NOTICE "************************************************************")
  message(NOTICE "================================")
  message(NOTICE "OOMPH-LIB CONFIGURATION OPTIONS:")
  message(NOTICE "================================")
  foreach(OPTION IN LISTS OOMPH_CONFIG_VARS)
    message(NOTICE "  ${MARKER} ${OPTION}: '${${OPTION}}'")
  endforeach()

  # Get the list compile definitions
  get_directory_property(
    OOMPH_COMPILE_DEFINITIONS DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                                        COMPILE_DEFINITIONS)

  # Sort the list alphabetically
  list(SORT OOMPH_COMPILE_DEFINITIONS)

  # Now print the compile definitions
  message(NOTICE "")
  message(NOTICE "===============================")
  message(NOTICE "OOMPH-LIB COMPILER DEFINITIONS:")
  message(NOTICE "===============================")
  foreach(DEFN ${OOMPH_COMPILE_DEFINITIONS})
    message(NOTICE "  ${MARKER} ${DEFN}")
  endforeach()
  message(NOTICE "************************************************************")
  message(NOTICE "")
endfunction()
# ------------------------------------------------------------------------------
