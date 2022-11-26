# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE
        "Entered linearised_free_surface_axisym_navier_stokes subdirectory")

set(SUBDIRS two_layer_interface_nonaxisym_perturbations)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE
        "Leaving linearised_free_surface_axisym_navier_stokes subdirectory")
# ------------------------------------------------------------------------------
