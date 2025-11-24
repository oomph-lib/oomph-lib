# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered linear_elasticity subdirectory")

set(SUBDIRS periodic_load refineable_periodic_load)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving linear_elasticity subdirectory")
# ------------------------------------------------------------------------------
