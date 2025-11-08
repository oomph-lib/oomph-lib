# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered pml_fourier_decomposed_helmholtz subdirectory")

set(SUBDIRS oscillating_sphere)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving pml_fourier_decomposed_helmholtz subdirectory")
# ------------------------------------------------------------------------------
