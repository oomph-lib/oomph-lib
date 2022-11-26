# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered biharmonic subdirectory")

set(SUBDIRS two_d_biharmonic)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving biharmonic subdirectory")
# ------------------------------------------------------------------------------
