# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered helmholtz subdirectory")

set(SUBDIRS scattering unstructured_scattering)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving helmholtz subdirectory")
# ------------------------------------------------------------------------------
