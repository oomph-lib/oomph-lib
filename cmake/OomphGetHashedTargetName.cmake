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
#     oomph_get_hashed_target_name(<EXECUTABLE_NAME> <OUTPUT_VARIABLE>)
#
# Example:
#
#     # Define executable
#     oomph_add_executable(NAME one_d_poisson
#                          SOURCES one_d_poisson.cc
#                          LIBRARIES oomph::poisson)
#
#     # Modify the target C++ standard
#     oomph_get_hashed_target_name(one_d_poisson HASHED_TARGET_NAME)
#     set_target_properties(${HASHED_TARGET_NAME} PROPERTIES CXX_STANDARD 20)
#
# It is worth noting that the argument supplied to CXX_DEFINITIONS does not
# require a -D prefix (as is usually required to indicate that it is a
# preprocessor definition), but you can supply it if you wish.
# =============================================================================
# cmake-format: on
include_guard()

function(oomph_get_hashed_target_name EXECUTABLE_NAME HASHED_TARGET_NAME)
  string(SHA1 PATH_HASH "${CMAKE_CURRENT_SOURCE_DIR}")
  string(SUBSTRING ${PATH_HASH} 0 7 PATH_HASH)
  set(${HASHED_TARGET_NAME} ${EXECUTABLE_NAME}_${PATH_HASH} PARENT_SCOPE)
endfunction()
# ------------------------------------------------------------------------------
