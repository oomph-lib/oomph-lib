# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered young_laplace subdirectory")

set(SUBDIRS young_laplace contact_angle)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving young_laplace subdirectory")
# ------------------------------------------------------------------------------
