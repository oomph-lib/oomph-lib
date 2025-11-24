# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered meshes subdirectory")

set(SUBDIRS
    mesh_list
    third_party_meshes
    mesh_from_xfig
    mesh_from_inline_triangle
    mesh_from_inline_triangle_internal_boundaries
    mesh_from_triangle
    mesh_from_vmtk
    mesh_from_geompack
    mesh_from_tetgen)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving meshes subdirectory")
# ------------------------------------------------------------------------------
