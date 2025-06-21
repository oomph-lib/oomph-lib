# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered fourier_decomposed_acoustic_fsi subdirectory")

set(SUBDIRS sphere unstructured_sphere)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving fourier_decomposed_acoustic_fsi subdirectory")
# ------------------------------------------------------------------------------
