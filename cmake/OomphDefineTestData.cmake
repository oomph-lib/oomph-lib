# cmake-format: off
# =============================================================================
# DESCRIPTION:
# ------------
# Defines the data required for tests in the current directory and the rule for
# copying it into the build directory.
#
# USAGE:
# ------
#   oomph_define_test_data(<file-or-dir-1> <file-or-dir-2> ... <file-or-dir-N>)
#
# EXAMPLE:
# --------
#   oomph_define_test_data(validata validate.sh my_extra_data_file.dat)
#
# NOTE: Arguments to DEPENDS_ON must be already-defined executables
# or targets (i.e. defined via add_executable() or oomph_add_executable()).
# =============================================================================
# cmake-format: on
include_guard()

# ------------------------------------------------------------------------------
function(oomph_define_test_data)
  # Hash the path to create a unique ID for our targets but shorten it to the
  # first 7 characters for brevity. A unique ID is required to avoid clashes
  # with targets in other directories
  string(SHA1 PATH_HASH "${CMAKE_CURRENT_SOURCE_DIR}")
  string(SUBSTRING ${PATH_HASH} 0 7 PATH_HASH)

  set(REQUIREMENTS_WITH_PATHS)
  set(TEST_BYPRODUCTS)

  # Add on the extra requirements
  foreach(REQUIREMENT IN LISTS ARGN)
    list(APPEND REQUIREMENTS_WITH_PATHS
         "${CMAKE_CURRENT_LIST_DIR}/${REQUIREMENT}")
    list(APPEND TEST_BYPRODUCTS "${CMAKE_CURRENT_BINARY_DIR}/${REQUIREMENT}")
  endforeach()

  # Declare a copy_... target for this directory if it hasn't been defined
  # already. We'll use this target to copy the required files to the build
  # directory
  if(NOT TARGET copy_${PATH_HASH})
    add_custom_target(copy_${PATH_HASH} WORKING_DIRECTORY "${CMAKE_BINARY_DIR}")
  endif()

  # Flag used to control whether files are symlinked instead of copied; keeping
  # the option to copy files around just in case we need it later on (but I
  # doubt it)
  set(SYMLINK_TEST_FILES_INSTEAD_OF_COPY TRUE)

  # Add each requirement to the copy target as a file-copy command or as a
  # directory-copy command. All of these commands will be executed when the
  # copy_<path-hash> target is called
  foreach(REQUIREMENT IN LISTS REQUIREMENTS_WITH_PATHS)
    if(SYMLINK_TEST_FILES_INSTEAD_OF_COPY)
      add_custom_command(
        TARGET copy_${PATH_HASH}
        POST_BUILD
        COMMAND ln -sf "${REQUIREMENT}" "${CMAKE_CURRENT_BINARY_DIR}")
    else()
      if(IS_DIRECTORY "${REQUIREMENT}")
        add_custom_command(
          TARGET copy_${PATH_HASH}
          POST_BUILD
          COMMAND cp -ur "${REQUIREMENT}" "${CMAKE_CURRENT_BINARY_DIR}")
      else()
        add_custom_command(
          TARGET copy_${PATH_HASH}
          POST_BUILD
          COMMAND ${CMAKE_COMMAND} -E copy_if_different "${REQUIREMENT}"
                  "${CMAKE_CURRENT_BINARY_DIR}")
      endif()
    endif()
  endforeach()

  # Identify the files that we'll copy as by-products so that they can be
  # cleaned up by running "make clean" if the user uses Makefile Generators
  add_custom_command(
    TARGET copy_${PATH_HASH}
    POST_BUILD
    BYPRODUCTS ${TEST_BYPRODUCTS})
endfunction()
# ------------------------------------------------------------------------------
