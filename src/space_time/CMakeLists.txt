# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered space_time subdirectory")

# The subdirectories to build
set(SRC_SUBDIRS space_time_block_preconditioner space_time_navier_stokes
    space_time_unsteady_heat)

# Add each subdirectory target to the build
foreach(SUBDIR ${SRC_SUBDIRS})
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving space_time subdirectory")
# ------------------------------------------------------------------------------
