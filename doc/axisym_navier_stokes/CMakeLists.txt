# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered axisym_navier_stokes subdirectory")

set(SUBDIRS spin_up axi_static_cap two_layer_interface_axisym)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving axisym_navier_stokes subdirectory")
# ------------------------------------------------------------------------------
