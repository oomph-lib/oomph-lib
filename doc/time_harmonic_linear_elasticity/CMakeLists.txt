# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered time_harmonic_linear_elasticity subdirectory")

set(SUBDIRS elastic_annulus unstructured_elastic_annulus)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving time_harmonic_linear_elasticity subdirectory")
# ------------------------------------------------------------------------------
