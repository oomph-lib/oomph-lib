# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE
        "Entered generalised_newtonian_axisym_navier_stokes subdirectory")

set(SUBDIRS axisym_vibrating_shell)

foreach(SUBDIR ${SUBDIRS})
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE
        "Leaving generalised_newtonian_axisym_navier_stokes subdirectory")
# ------------------------------------------------------------------------------
