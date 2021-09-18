# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered axisym_advection_diffusion subdirectory")

set(SUBDIRS pipe)

foreach(SUBDIR ${SUBDIRS})
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving axisym_advection_diffusion subdirectory")
# ------------------------------------------------------------------------------
