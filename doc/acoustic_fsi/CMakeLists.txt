# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered acoustic_fsi subdirectory")

set(SUBDIRS acoustic_fsi_annulus unstructured_acoustic_fsi_annulus)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving acoustic_fsi subdirectory")
# ------------------------------------------------------------------------------
