# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered preconditioners subdirectory")

set(SUBDIRS lsc_navier_stokes lgr_navier_stokes fsi prescribed_displ_lagr_mult
    pseudo_solid_fsi)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving preconditioners subdirectory")
# ------------------------------------------------------------------------------
