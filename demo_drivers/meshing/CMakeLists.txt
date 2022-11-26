# ------------------------------------------------------------------------------
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered meshing subdirectory")

set(SUBDIRS
    mesh_from_inline_triangle
    mesh_from_inline_triangle_internal_boundaries
    mesh_from_triangle
    mesh_from_geompack
    mesh_from_tetgen
    mesh_from_xfig_triangle
    mesh_from_vmtk
    quad_from_triangle_mesh
    quadratic_vmtk_tetgen
    adaptive_tet_meshes)

foreach(SUBDIR IN LISTS SUBDIRS)
  add_subdirectory(${SUBDIR})
endforeach()

message(VERBOSE "Leaving meshing subdirectory")
# ------------------------------------------------------------------------------
