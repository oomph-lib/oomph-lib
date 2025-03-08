# cmake-format: off
# =============================================================================
# DESCRIPTION:
# ------------
# This script adds clang-tidy support to the library build process.
# When the flag OOMPH_ENABLE_CLANG_TIDY is enabled, this script:
#   - Searches for the clang-tidy executable.
#   - Fails the configuration with an error if clang-tidy is not found.
#   - Sets the CMAKE_CXX_CLANG_TIDY property to automatically run
#     clang-tidy on your C++ targets.
#
# USAGE:
# ------
# In your top-level CMakeLists.txt, add:
#
#   option(OOMPH_ENABLE_CLANG_TIDY "Enable clang-tidy support for static analysis" OFF)
#   include(OomphAddClangTidySupport)
#
# assuming this cmake/ folder is in your CMAKE_MODULE_PATH.
#
# When OOMPH_ENABLE_CLANG_TIDY is ON, clang-tidy will be used for static analysis.
#
# -----------------------------------------------------------------------------
# IMPORTANT:
# ----------
# Make sure that clang-tidy is installed and available in your system's PATH.
# If clang-tidy cannot be found when OOMPH_ENABLE_CLANG_TIDY is enabled,
# the configuration process will terminate with a fatal error.
# =============================================================================
include_guard()
# cmake-format: on
# ------------------------------------------------------------------------------
# Check if clang-tidy support is enabled.
if(OOMPH_ENABLE_CLANG_TIDY)
  if(NOT EXISTS "${CMAKE_SOURCE_DIR}/.clang-tidy")
    message(
      STATUS
        "${CMAKE_SOURCE_DIR}/.clang-tidy file not found. Skipping clang-tidy integration."
    )
  else()
    message(STATUS "clang-tidy support enabled. Searching for clang-tidy...")

    # Find the clang-tidy executable.
    find_program(CLANG_TIDY_EXE clang-tidy REQUIRED)

    # Tell the user where we found it
    message(STATUS "Found clang-tidy: ${CLANG_TIDY_EXE}")
    message(
      STATUS "Found .clang-tidy file file: ${CMAKE_SOURCE_DIR}/.clang-tidy")

    # Optionally define additional clang-tidy options (e.g., specific checks).
    # Modify CLANG_TIDY_OPTIONS as needed.
    set(OOMPH_CLANG_TIDY_OPTIONS "")

    # Set clang-tidy for all C++ targets via the CMAKE_CXX_CLANG_TIDY property.
    set(CMAKE_CXX_CLANG_TIDY "${CLANG_TIDY_EXE}")
  endif()
endif()
# ------------------------------------------------------------------------------
