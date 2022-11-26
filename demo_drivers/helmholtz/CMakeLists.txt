# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered helmholtz subdirectory")

set(SUBDIRS scattering point_source)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving helmholtz subdirectory")
# ------------------------------------------------------------------------------
