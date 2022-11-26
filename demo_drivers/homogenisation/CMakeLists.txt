# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered homogenisation subdirectory")

set(SUBDIRS two_dim)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving homogenisation subdirectory")
# ------------------------------------------------------------------------------
