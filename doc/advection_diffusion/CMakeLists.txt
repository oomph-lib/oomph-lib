# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered advection_diffusion subdirectory")

set(SUBDIRS two_d_adv_diff_adapt two_d_adv_diff_flux_bc two_d_adv_diff_SUPG)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving advection_diffusion subdirectory")
# ------------------------------------------------------------------------------
