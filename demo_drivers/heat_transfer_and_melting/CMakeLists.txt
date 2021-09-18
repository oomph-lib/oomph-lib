# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered heat_transfer_and_melting subdirectory")

set(SUBDIRS two_d_unsteady_heat_melt)

foreach(SUBDIR ${SUBDIRS})
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving heat_transfer_and_melting subdirectory")
# ------------------------------------------------------------------------------
