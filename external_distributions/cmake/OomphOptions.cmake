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
#     FLAG OOMPH_THIRD_PARTY_INSTALL_DIR
#     DEFAULT "${CMAKE_CURRENT_LIST_DIR}/install"
#     DOCSTRING "Base installation directory for third-party libraries."
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
  set(FLAGS IS_REQUIRED OVERWRITE)
  set(SINGLE_VALUE_ARGS FLAG DEFAULT DOCSTRING)
  set(MULTI_VALUE_ARGS)

  # Process the arguments passed in
  cmake_parse_arguments(PARSE_ARGV 0 ${PREFIX} "${FLAGS}"
                        "${SINGLE_VALUE_ARGS}" "${MULTI_VALUE_ARGS}")

  # Redefine the variables in this scope without a prefix for clarity
  set(FLAG ${${PREFIX}_FLAG})
  set(DEFAULT ${${PREFIX}_DEFAULT})
  set(DOCSTRING ${${PREFIX}_DOCSTRING})
  set(IS_REQUIRED ${${PREFIX}_IS_REQUIRED})
  set(OVERWRITE ${${PREFIX}_OVERWRITE})

  # Would be weird to provide a default argument but require an input
  if(DEFAULT AND IS_REQUIRED)
    message(
      FATAL_ERROR
        "Can't provide oomph_path_option(...) a default value and also set 'IS_REQUIRED'!"
    )
  endif()

  # Forceably overwrite the cache variable everytime the configure step is run?
  set(EXTRA_FLAGS)
  if(OVERWRITE)
    set(EXTRA_FLAGS FORCE)
  endif()

  # Create a cache
  if(${FLAG})
    set(${FLAG} ${${FLAG}} CACHE PATH "${DOCSTRING}")
  elseif(DEFAULT)
    set(${FLAG} ${DEFAULT} CACHE PATH "${DOCSTRING}")
  elseif(IS_REQUIRED)
    message(
      FATAL_ERROR
        "Argument '${FLAG}' is required but you did not specify a value! Set it using -D${FLAG}=\"...\" at the commandline."
    )
  endif()

  # Update the list of commandline flags that we take but remember to wipe
  # duplicates so we don't keep appending the same variable again and again
  set(OOMPH_CONFIG_VARS ${OOMPH_CONFIG_VARS} ${FLAG})
  list(REMOVE_DUPLICATES OOMPH_CONFIG_VARS)
  set(OOMPH_CONFIG_VARS ${OOMPH_CONFIG_VARS} CACHE INTERNAL
      "List of oomph-lib configuration options" FORCE)
endfunction()
# ------------------------------------------------------------------------------
