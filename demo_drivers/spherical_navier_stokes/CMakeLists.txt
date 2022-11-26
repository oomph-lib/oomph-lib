# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered spherical_navier_stokes subdirectory")

set(SUBDIRS spin_up spherical_couette steady_rot)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving spherical_navier_stokes subdirectory")
# ------------------------------------------------------------------------------
