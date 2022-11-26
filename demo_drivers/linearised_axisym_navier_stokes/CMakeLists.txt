# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered linearised_axisym_navier_stokes subdirectory")

set(SUBDIRS counter_rotating_disks time_periodic_taylor_couette)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving linearised_axisym_navier_stokes subdirectory")
# ------------------------------------------------------------------------------
