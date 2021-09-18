# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered polar_navier_stokes subdirectory")

set(SUBDIRS jeffery_hamel)

foreach(SUBDIR ${SUBDIRS})
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving polar_navier_stokes subdirectory")
# ------------------------------------------------------------------------------
