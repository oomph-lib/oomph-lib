# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered axisym_linear_elasticity subdirectory")

set(SUBDIRS cylinder)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving axisym_linear_elasticity subdirectory")
# ------------------------------------------------------------------------------
