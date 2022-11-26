# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered disable_ALE subdirectory")

set(SUBDIRS unsteady_heat navier_stokes)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving disable_ALE subdirectory")
# ------------------------------------------------------------------------------
