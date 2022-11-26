# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered womersley subdirectory")

set(SUBDIRS one_d_womersley two_d_womersley)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving womersley subdirectory")
# ------------------------------------------------------------------------------
