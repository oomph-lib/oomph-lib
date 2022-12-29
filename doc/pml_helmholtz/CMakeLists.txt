# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered pml_helmholtz subdirectory")

set(SUBDIRS scattering)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving pml_helmholtz subdirectory")
# ------------------------------------------------------------------------------
