# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered axisym_poroelasticity subdirectory")

set(SUBDIRS unstructured_two_d_curved)

foreach(SUBDIR ${SUBDIRS})
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving axisym_poroelasticity subdirectory")
# ------------------------------------------------------------------------------
