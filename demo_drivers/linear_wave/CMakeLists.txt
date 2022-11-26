# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered linear_wave subdirectory")

set(SUBDIRS two_d_linear_wave two_d_linear_wave_adapt)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving linear_wave subdirectory")
# ------------------------------------------------------------------------------
