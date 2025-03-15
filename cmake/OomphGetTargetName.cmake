# cmake-format: off
# =============================================================================
# DESCRIPTION:
# ------------
# Get the target name associated with the specified executable.
#
# oomph-lib has a lot of targets with similar names. To avoid clashes when
# building, e.g. all of the demo drivers, we append a 7-character hash to
# the executable name to construct the target name. (Note that we make sure
# to keep the name of the output executable the same as the one that is
# specified when calling, e.g. oomph_add_executable.)
#
# USAGE:
# ------
#     oomph_get_target_name(<EXECUTABLE_NAME> <OUTPUT_VARIABLE>)
#
# EXAMPLE:
# --------
#
#     # Define executable
#     oomph_add_executable(NAME one_d_poisson
#                          SOURCES one_d_poisson.cc
#                          LIBRARIES oomph::poisson)
#
#     # Modify the target C++ standard
#     oomph_get_target_name(one_d_poisson OOMPH_TARGET_NAME)
#     set_target_properties(${OOMPH_TARGET_NAME} PROPERTIES CXX_STANDARD 20)
#
# =============================================================================
include_guard()

# ------------------------------------------------------------------------------
function(oomph_get_target_name EXEC_NAME OOMPH_TARGET_NAME)
  string(SHA1 PATH_HASH "${CMAKE_CURRENT_LIST_DIR}")
  string(SUBSTRING ${PATH_HASH} 0 7 PATH_HASH)
  set(${OOMPH_TARGET_NAME} ${EXEC_NAME}_${PATH_HASH} PARENT_SCOPE)
endfunction()
# ------------------------------------------------------------------------------
# cmake-format: on
