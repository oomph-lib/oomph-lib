# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered shell subdirectory")

set(SUBDIRS clamped_shell clamped_or_pinned_shell oscillating_shell plate)

foreach(SUBDIR ${SUBDIRS})
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving shell subdirectory")
# ------------------------------------------------------------------------------
