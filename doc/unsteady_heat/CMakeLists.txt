# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered unsteady_heat subdirectory")

set(SUBDIRS
    two_d_unsteady_heat
    two_d_unsteady_heat2
    two_d_unsteady_heat_t_adapt
    two_d_unsteady_heat_adapt
    two_d_unsteady_heat_ALE
    two_d_unsteady_heat_2adapt)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving unsteady_heat subdirectory")
# ------------------------------------------------------------------------------
