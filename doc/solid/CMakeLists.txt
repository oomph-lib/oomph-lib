# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered solid subdirectory")

set(SUBDIRS
    static_fish
    disk_compression
    disk_oscillation
    shock_disk
    airy_cantilever
    solid_theory
    prescribed_displ_lagr_mult
    three_d_cantilever
    simple_shear
    vmtk_solid
    unstructured_solid
    unstructured_three_d_solid
    unstructured_adaptive_solid
    compressed_square)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving solid subdirectory")
# ------------------------------------------------------------------------------
