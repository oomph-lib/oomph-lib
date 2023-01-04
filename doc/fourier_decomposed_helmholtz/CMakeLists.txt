# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered fourier_decomposed_helmholtz subdirectory")

set(SUBDIRS sphere_scattering adaptive_sphere_scattering)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving fourier_decomposed_helmholtz subdirectory")
# ------------------------------------------------------------------------------
