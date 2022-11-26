# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered eigenproblems subdirectory")

set(SUBDIRS)

if(HAVE_ARPACK_SOURCES)
  list(APPEND SUBDIRS harmonic orr_sommerfeld)
endif()

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving eigenproblems subdirectory")
# ------------------------------------------------------------------------------
