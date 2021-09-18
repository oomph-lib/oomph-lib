# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered multi_domain subdirectory")

set(SUBDIRS
    boussinesq_convection
    fsi_collapsible_channel
    fsi_osc_ring
    fsi_channel_with_leaflet
    turek_flag
    pseudo_solid_collapsible_tube)

foreach(SUBDIR ${SUBDIRS})
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving multi_domain subdirectory")
# ------------------------------------------------------------------------------
