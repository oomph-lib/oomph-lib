# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered gzip subdirectory")

set(SUBDIRS one_d_poisson)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving gzip subdirectory")
# ------------------------------------------------------------------------------
