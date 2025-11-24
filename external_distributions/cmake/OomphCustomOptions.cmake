# cmake-format: off
# =============================================================================
# DESCRIPTION:
# ------------
#
# USAGE:
# ------
#
#   oomph_option(OOMPH_ENABLE_MPI "Enable the use of MPI for parallel processing" OFF)
#
#   oomph_path_option(
#     FLAG OOMPH_USE_OPENBLAS_FROM
#     DOCSTRING "The path to a preinstalled version of OpenBLAS."
#     REQUIRED
#   )
#
# =============================================================================
# cmake-format: on
include_guard()
include(CMakeParseArguments)

# ------------------------------------------------------------------------------
# cmake-format: off
function(oomph_option VAR DOCSTRING VALUE)
  option(${VAR} "${DOCSTRING}" ${VALUE})
  set(OOMPH_CONFIG_VARS ${OOMPH_CONFIG_VARS} ${VAR} CACHE INTERNAL
      "List of oomph-lib configuration options" FORCE)
endfunction()
# cmake-format: on
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
function(oomph_path_option)
  set(PREFIX ARG)
  set(FLAGS REQUIRED)
  set(SINGLE_VALUE_ARGS FLAG DEFAULT DOCSTRING)
  set(MULTI_VALUE_ARGS)

  # Process the arguments passed in
  cmake_parse_arguments(PARSE_ARGV 0 ${PREFIX} "${FLAGS}"
                        "${SINGLE_VALUE_ARGS}" "${MULTI_VALUE_ARGS}")

  # Redefine the variables in this scope without a prefix for clarity
  set(FLAG ${${PREFIX}_FLAG})
  set(DEFAULT "${${PREFIX}_DEFAULT}") # defaults to empty string
  set(DOCSTRING "${${PREFIX}_DOCSTRING}")
  set(REQUIRED ${${PREFIX}_REQUIRED})

  if(NOT FLAG)
    message(FATAL_ERROR "oomph_path_option() requires FLAG <varname>.")
  endif()

  # Create a cache entry only if none exists yet (so it shows in ccmake)
  if(NOT DEFINED CACHE{${FLAG}})
    if("${${FLAG}}" STREQUAL "")
      set(${FLAG} "${DEFAULT}" CACHE PATH "${DOCSTRING}")
    else()
      set(${FLAG} "${${FLAG}}" CACHE PATH "${DOCSTRING}")
    endif()
  endif()

  # If the path variable is required, make sure it's a non-empty string
  if(REQUIRED AND ("${${FLAG}}" STREQUAL ""))
    message(
      FATAL_ERROR
        "Argument '${FLAG}' is required but you did not specify a non-empty value! "
        "Set it using -D${FLAG}=\"...\" at the commandline.")
  endif()

  # Track config vars (deduped)
  set(OOMPH_CONFIG_VARS ${OOMPH_CONFIG_VARS} ${FLAG})
  list(REMOVE_DUPLICATES OOMPH_CONFIG_VARS)
  set(OOMPH_CONFIG_VARS ${OOMPH_CONFIG_VARS} CACHE INTERNAL
      "List of oomph-lib configuration options" FORCE)
endfunction()
# ------------------------------------------------------------------------------
