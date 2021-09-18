# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered multigrid subdirectory")

set(SUBDIRS scalar_multigrid helmholtz_multigrid)

foreach(SUBDIR ${SUBDIRS})
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving multigrid subdirectory")
# ------------------------------------------------------------------------------
