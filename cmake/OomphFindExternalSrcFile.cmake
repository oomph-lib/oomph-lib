# cmake-format: off
# =============================================================================
# DESCRIPTION:
# ------------
# Find a file FILE required by external_src. We search for the file in current
# directory first, then in the paths specified by SEARCH_PATHS (if provided).
# If the file is found in one of the search paths then a symlink will be created
# to that file. Otherwise, we use the file DUMMY_FILE (if provided) in its
# place. The chosen file will then be appended to the TARGET argument. If the
# file wasn't found, a dummy file wasn't provided and the REQUIRED flag was set,
# then we raise a FATAL_ERROR. Otherwise, we do nothing and continue quietly.
#
# USAGE:
# ------
# include(OomphFindExternalSrcFile)
# oomph_find_external_src_file(FILE          <file-to-locate-or-symlink>
#                              [TARGET       <variable-to-append-file-to>]
#                              [SEARCH_PATHS <additional-paths-to-search>]
#                              [SET_FLAG     <boolean-variable-to-set>]
#                              [DUMMY_FILE   <headers-to-be-symlinked-to>]
#                              [REQUIRED]
#                              [TOUCH_DUMMY_FILE])
#
# The variable SET_FLAG will be set and made visible across the entire project.
# As an example, this may be set to HAVE_ARPACK_SOURCES. If the REQUIRED flag
# specified then an error will be issued if the file couldn't be found and no
# dummy file was provided. If the TOUCH_DUMMY_FILE flag is specified then the
# dummy file will be "touch"-ed (i.e. it has the same effect as the "touch"
# shell command).
# =============================================================================
# cmake-format: on
include_guard()

function(oomph_find_external_src_file)
  # Define the supported set of keywords. The value of PREFIX will be prepended
  # to each keyword after passing through cmake_parse_arguments(...).
  set(PREFIX ARG)
  set(FLAGS REQUIRED TOUCH_DUMMY_FILE)
  set(SINGLE_VALUE_ARGS FILE TARGET SET_FLAG DUMMY_FILE)
  set(MULTI_VALUE_ARGS SEARCH_PATHS)

  # Process the arguments passed in
  include(CMakeParseArguments)
  cmake_parse_arguments(PARSE_ARGV 0 ${PREFIX} "${FLAGS}"
                        "${SINGLE_VALUE_ARGS}" "${MULTI_VALUE_ARGS}")

  # Store the parsed keyword arguments in variables with clearer names
  set(FILE ${${PREFIX}_FILE})
  set(TARGET ${${PREFIX}_TARGET})
  set(REQUIRED ${${PREFIX}_REQUIRED})
  set(SET_FLAG ${${PREFIX}_SET_FLAG})
  set(DUMMY_FILE ${${PREFIX}_DUMMY_FILE})
  set(SEARCH_PATHS ${${PREFIX}_SEARCH_PATHS})
  set(TOUCH_DUMMY_FILE ${${PREFIX}_TOUCH_DUMMY_FILE})

  # Assume the file hasn't been found by default and set the flag to false/no
  if(SET_FLAG)
    set(${SET_FLAG} NO)
  endif()

  # Search for the file in the current directory and in the private external
  # sources directory. Don't look anywhere else (NO_DEFAULT_PATH)
  find_file(
    PATH_TO_FILE ${FILE}
    PATHS "${CMAKE_CURRENT_LIST_DIR}" "${SEARCH_PATHS}"
    NO_DEFAULT_PATH)

  # If we were able to find the file
  if(PATH_TO_FILE)
    # Tell the user
    file(RELATIVE_PATH REL_PATH ${CMAKE_SOURCE_DIR} "${PATH_TO_FILE}")
    message(VERBOSE "Found ${REL_PATH}")

    # Mark the file as being found
    if(SET_FLAG)
      set(${SET_FLAG} YES)
    endif()

    # If the file isn't in the calling directory, symlink it
    cmake_path(GET PATH_TO_FILE PARENT_PATH FILE_LOCATION)
    if(NOT FILE_LOCATION STREQUAL ${CMAKE_CURRENT_LIST_DIR})
      message(VERBOSE "Creating symlink.")
      file(CREATE_LINK ${FILE_LOCATION} "${CMAKE_CURRENT_LIST_DIR}/${FILE}"
           SYMBOLIC)
    endif()
  elseif(DUMMY_FILE)
    # If the file couldn't be found use the dummy version (if provided)
    set(OUTPUT_STR "Couldn't find ${FILE}. ")

    # Create an empty dummy file if the file doesn't exist already
    if(TOUCH_DUMMY_FILE)
      string(APPEND OUTPUT_STR "Creating ${DUMMY_FILE} dummy file instead.")
      file(TOUCH ${DUMMY_FILE})
    else()
      string(APPEND OUTPUT_STR "Using ${DUMMY_FILE} instead.")
    endif()

    message(VERBOSE "${OUTPUT_STR}")
    set(PATH_TO_FILE ${DUMMY_FILE})
  else()
    # The file couldn't be found. Issue an error if the file was a requirement
    if(REQUIRED)
      message(FATAL_ERROR "Couldn't find ${FILE}!")
    else()
      message(STATUS "Couldn't find ${FILE}. Doing nothing.")
    endif()

    # Finish here
    return()
  endif()

  # If the user wanted the filename appended to a specified target
  if(TARGET)
    # Append the filename of the chosen file
    cmake_path(GET PATH_TO_FILE FILENAME FILENAME)
    set(UPDATED_TARGET ${${TARGET}} ${FILENAME})

    # Update the specified target and make it visible in the parent scope
    set(${TARGET} ${UPDATED_TARGET} PARENT_SCOPE)
  endif()

  # If the user wanted a flag to be set and made visible to the entire project
  # scope then set it here. For this we use a cache variable.
  if(SET_FLAG)
    set(${SET_FLAG} ${${SET_FLAG}} CACHE INTERNAL "")
  endif()

  # Clean up: delete cache variables that need to be set each time the function
  # is called. If these variables aren't cleaned up then the result from the
  # previous call will be used (which we don't want!)
  unset(PATH_TO_FILE CACHE)
endfunction()
