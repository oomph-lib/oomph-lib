# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered eigenproblems subdirectory")

set(SUBDIRS harmonic complex_harmonic)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving eigenproblems subdirectory")
# ------------------------------------------------------------------------------
