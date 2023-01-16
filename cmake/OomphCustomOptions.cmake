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
# cmake-format: on
