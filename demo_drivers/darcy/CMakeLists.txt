# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered darcy subdirectory")

set(SUBDIRS unstructured_two_d_circle)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving darcy subdirectory")
# ------------------------------------------------------------------------------
