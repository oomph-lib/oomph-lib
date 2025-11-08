# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered multi_physics subdirectory")

set(SUBDIRS b_convection refine_b_convect multi_domain_ref_b_convect thermo
    rayleigh_instability_surfactant)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving multi_physics subdirectory")
# ------------------------------------------------------------------------------
