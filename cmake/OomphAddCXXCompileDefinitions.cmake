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

function(oomph_add_cxx_compile_definitions)
  foreach(CXX_DEFINITION IN LISTS ARGN)
    add_compile_definitions($<$<COMPILE_LANGUAGE:CXX>:${CXX_DEFINITION}>)
  endforeach()
endfunction()

function(oomph_add_c_compile_definitions)
  foreach(C_DEFINITION IN LISTS ARGN)
    add_compile_definitions($<$<COMPILE_LANGUAGE:C>:${C_DEFINITION}>)
  endforeach()
endfunction()
# ------------------------------------------------------------------------------
