# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered beam subdirectory")

set(SUBDIRS tensioned_string steady_ring steady_third_ring unsteady_ring)

foreach(SUBDIR ${SUBDIRS})
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving beam subdirectory")
# ------------------------------------------------------------------------------
