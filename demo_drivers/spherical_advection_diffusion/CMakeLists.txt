# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered spherical_advection_diffusion subdirectory")

set(SUBDIRS eluting_sphere)

foreach(SUBDIR ${SUBDIRS})
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving spherical_advection_diffusion subdirectory")
# ------------------------------------------------------------------------------
