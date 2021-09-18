# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(
  VERBOSE
  "Entered time_harmonic_fourier_decomposed_linear_elasticity subdirectory")

set(SUBDIRS cylinder)

foreach(SUBDIR ${SUBDIRS})
  add_subdirectory(${SUBDIR})
endforeach()

message(
  VERBOSE
  "Leaving time_harmonic_fourier_decomposed_linear_elasticity subdirectory")
# ------------------------------------------------------------------------------
