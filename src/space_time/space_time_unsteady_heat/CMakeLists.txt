# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered space_time_unsteady_heat subdirectory")

# The subdirectories to build
set(SRC_SUBDIRS galerkin_equal_order galerkin_petrov_equal_order
    galerkin_petrov_mixed_order)

# Add each subdirectory target to the build
foreach(SUBDIR ${SRC_SUBDIRS})
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving space_time_unsteady_heat subdirectory")
# ------------------------------------------------------------------------------
