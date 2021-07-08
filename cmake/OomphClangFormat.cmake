# cmake-format: off
# =============================================================================
# DESCRIPTION:
# ------------
# Provides automatic formatting of C++ code. Need to confine its use to
# oomph-lib sources only.
#
# NOT FINISHED OFF.
#
# USAGE:
# ------
#     include(OomphClangFormat)
#
# =============================================================================
# cmake-format: on

if(OOMPH_ENABLE_CLANG_FORMAT)
  # Additional targets to perform clang-format/clang-tidy; get all project files
  file(GLOB_RECURSE ALL_CXX_SOURCE_FILES *.[cht]pp *.cc *.h)

  # Search the system paths for the clang-tidy tool
  find_program(CLANG_FORMAT_PROGRAM "clang-format")

  # Run clang-format during the build if found
  if(CLANG_FORMAT_PROGRAM)
    add_custom_target(clang-format COMMAND ${CLANG_FORMAT_PROGRAM} -i
                                           -style=file ${ALL_CXX_SOURCE_FILES})
  endif()
endif()
