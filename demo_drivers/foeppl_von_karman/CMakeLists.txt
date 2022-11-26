# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered foeppl_von_karman subdirectory")

set(SUBDIRS circular_disk)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving foeppl_von_karman subdirectory")
# ------------------------------------------------------------------------------
