# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered mpi subdirectory")

set(SUBDIRS distribution multi_domain solvers)

foreach(SUBDIR ${SUBDIRS})
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving mpi subdirectory")
# ------------------------------------------------------------------------------
