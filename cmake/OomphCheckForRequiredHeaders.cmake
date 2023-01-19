# cmake-format: off
# =============================================================================
# DESCRIPTION:
# ------------
# Make sure we have certain headers required for the oomph-lib build.
#
# USAGE:
# ------
# include(OomphCheckForRequiredHeaders)
# =============================================================================
# cmake-format: on
include_guard()

# Use built-in CMake module to search for headers
include(CheckIncludeFiles)

# The required headers
set(REQUIRED_HEADERS limits.h stdlib.h string.h sys/time.h)

# Use the CheckIncludeFiles module to see if we can find must-have headers
check_include_files("${REQUIRED_HEADERS}" OOMPH_HAS_REQUIRED_HEADERS)

# Issue an error if the headers weren't found
if(NOT OOMPH_HAS_REQUIRED_HEADERS)
  message(FATAL_ERROR "Trouble -- unable to find required header(s)!")
endif()

# Make sure we have the desired headers (i.e. not mandatory)
check_include_files(execinfo.h OOMPH_HAS_STACKTRACE)
check_include_files(unistd.h OOMPH_HAS_UNISTDH)
check_include_files(malloc.h OOMPH_HAS_MALLOCH)

# Add compiler flags
if(OOMPH_HAS_STACKTRACE)
  oomph_add_cxx_compile_definitions(OOMPH_HAS_STACKTRACE)
endif()
if(OOMPH_HAS_UNISTDH)
  oomph_add_cxx_compile_definitions(OOMPH_HAS_UNISTDH)
endif()
if(OOMPH_HAS_MALLOCH)
  oomph_add_cxx_compile_definitions(OOMPH_HAS_MALLOCH)
endif()
