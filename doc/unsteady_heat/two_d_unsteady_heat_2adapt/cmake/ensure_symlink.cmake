# ------------------------------------------------------------------------------
# Usage:  cmake -P ensure_symlink.cmake <SOURCE> <DEST>
# ------------------------------------------------------------------------------
if(NOT "${CMAKE_ARGC}" EQUAL 5)
  message(FATAL_ERROR "Usage: ensure_symlink.cmake <SOURCE> <DEST>")
endif()

math(EXPR SRC_IDX "${CMAKE_ARGC} - 2") # index of second-to-last argument
math(EXPR DST_IDX "${CMAKE_ARGC} - 1") # index of last argument

set(SRC "${CMAKE_ARGV${SRC_IDX}}")
set(DST "${CMAKE_ARGV${DST_IDX}}")

# Compute the directory that will contain the symlink
get_filename_component(DST_DIR "${DST}" DIRECTORY)

# Compute the path from DST_DIR to SRC
file(RELATIVE_PATH REL_SRC "${DST_DIR}" "${SRC}")

if(NOT EXISTS "${DST}")
  # No file yet so create a link
  execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink "${REL_SRC}"
                          "${DST}")
elseif(NOT IS_SYMLINK "${DST}")
  # Exists but is a regular file/dir so remove + link
  file(REMOVE_RECURSE "${DST}")
  execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink "${REL_SRC}"
                          "${DST}")
else()
  file(READ_SYMLINK "${DST}" CURRENT_TARGET)
  if(NOT CURRENT_TARGET STREQUAL "${REL_SRC}")
    # Points somewhere else so relink
    file(REMOVE "${DST}")
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink "${REL_SRC}"
                            "${DST}")
  endif()
endif()
# ------------------------------------------------------------------------------
