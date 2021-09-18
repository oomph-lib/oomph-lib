# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered optimisation subdirectory")

set(SUBDIRS
    C_style_output
    disable_ALE
    linear_vs_nonlinear
    stored_shape_fcts
    use_coarse_base_meshes
    sparse_assemble
    fsi_jacobian_approximation)

foreach(SUBDIR ${SUBDIRS})
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving optimisation subdirectory")
# ------------------------------------------------------------------------------
