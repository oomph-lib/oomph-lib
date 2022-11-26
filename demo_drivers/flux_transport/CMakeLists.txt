# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered flux_transport subdirectory")

set(SUBDIRS advection euler)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving flux_transport subdirectory")
# ------------------------------------------------------------------------------
