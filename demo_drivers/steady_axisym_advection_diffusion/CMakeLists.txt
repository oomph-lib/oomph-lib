# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered steady_axisym_advection_diffusion subdirectory")

set(SUBDIRS steady_axisym_advection_diffusion)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving steady_axisym_advection_diffusion subdirectory")
# ------------------------------------------------------------------------------
