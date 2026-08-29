# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered axisym_solid subdirectory")

set(SUBDIRS spherical_shell_compression spherical_shell_unsteady)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving axisym_solid subdirectory")
# ------------------------------------------------------------------------------
