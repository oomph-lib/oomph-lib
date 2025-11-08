# cmake-format: off
# =============================================================================
# DESCRIPTION:
# ------------
# Helper to define compiler definitions for C++ files only.
#
# USAGE:
# ------
#     oomph_add_cxx_compile_definitions(<CXX-DEFINITION>)
#     oomph_add_cxx_compile_definitions(${<LIST-OF-CXX-DEFINITIONS>})
#
# EXAMPLE:
# --------
#     # Add -DPARANOID compiler definition when CMake is built in Debug mode
#     oomph_add_cxx_compile_definitions($<$<CONFIG:DEBUG,>:PARANOID>)
#
#     # Add -DOOMPH_HAS_MPI compiler definition
#     oomph_add_cxx_compile_definitions(OOMPH_HAS_MPI)
# =============================================================================
# cmake-format: on
include_guard()

# ------------------------------------------------------------------------------
function(oomph_add_cxx_compile_definitions)
  # Join all ARGN entries into a single string (in case the user quoted it)
  string(REPLACE ";" " " JOINED "${ARGN}")

  # Split on spaces into a list
  separate_arguments(DEF_LIST UNIX_COMMAND "${JOINED}")
  foreach(CXX_DEFINITION IN LISTS DEF_LIST)
    # Strip leading -D, if found
    string(REGEX REPLACE "^-D" "" CLEAN_DEF "${CXX_DEFINITION}")

    # Apply compile definition
    add_compile_definitions($<$<COMPILE_LANGUAGE:CXX>:${CLEAN_DEF}>)
  endforeach()
endfunction()

function(oomph_add_c_compile_definitions)
  # Join all ARGN entries into a single string (in case the user quoted it)
  string(REPLACE ";" " " JOINED "${ARGN}")

  # Split on spaces into a list
  separate_arguments(DEF_LIST UNIX_COMMAND "${JOINED}")
  foreach(C_DEFINITION IN LISTS DEF_LIST)
    # Strip leading -D, if found
    string(REGEX REPLACE "^-D" "" CLEAN_DEF "${C_DEFINITION}")

    # Apply compile definition
    add_compile_definitions($<$<COMPILE_LANGUAGE:C>:${CLEAN_DEF}>)
  endforeach()
endfunction()
# ------------------------------------------------------------------------------
