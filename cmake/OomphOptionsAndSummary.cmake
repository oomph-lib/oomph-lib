# cmake-format: off
# =============================================================================
# DESCRIPTION:
# ------------
#
# USAGE:
# ------
#
# =============================================================================
include_guard()

# ------------------------------------------------------------------------------
function(oomph_option VAR DOCSTRING VALUE)
  option(${VAR} "${DOCSTRING}" ${VALUE})
  set(OOMPH_CONFIG_VARS ${OOMPH_CONFIG_VARS} ${VAR} CACHE INTERNAL
    "List of oomph-lib configuration options" FORCE)
endfunction()
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
function(oomph_path_option VAR DOCSTRING VALUE)
  set(${VAR} ${VALUE} CACHE PATH "${DOCSTRING}")
  set(OOMPH_CONFIG_VARS ${OOMPH_CONFIG_VARS} ${VAR} CACHE INTERNAL
    "List of oomph-lib configuration options" FORCE)
endfunction()
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
function(oomph_print_config_summary)
  include(FeatureSummary)
  message("")
  message("************************************************************")
  message("OOMPH-LIB SETTINGS:")
  message("============================================================")
  foreach(OPTION IN LISTS OOMPH_CONFIG_VARS)
    message("  * ${OPTION}: '${${OPTION}}'")
  endforeach()
  message("************************************************************")
  message("")
endfunction()
# ------------------------------------------------------------------------------
# cmake-format: on
