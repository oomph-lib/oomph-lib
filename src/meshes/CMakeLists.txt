# ==============================================================================
# NOTE:
#
# Mesh files generally come in two parts: *.h and *.template.cc. Usually there
# are no *.cc files because meshes are fully templated. It therefore doesn't
# make any sense to compile mesh sources into a library. Instead, we combine the
# two files into a combined mesh header, *.h, which can be included into a
# driver code. This means that the entire code (headers and function
# implementations) can get (re-)compiled with the appropriate element types in
# each specific driver code. If the constant re-compilation of the (possibly
# lengthy) mesh sources in *.template.cc gets too costly, then *.h and
# *.template.cc can be included in a separate file where the instantiation of a
# mesh with a specific element type can be enforced with a template statement.
# ------------------------------------------------------------------------------
# Provide Domain definitions in their respective headers; don't write the
# definitions in source files. If you provide the definitions outside of the
# class, make sure to mark it as "inline" otherwise you may well get warnings
# about duplicate symbols!
# ==============================================================================
list(APPEND CMAKE_MESSAGE_INDENT " ")
message(VERBOSE "Entered meshes subdirectory")

# Define the templated headers
set(TEMPLATE_HEADERS
    one_d_mesh.h
    simple_rectangular_quadmesh.h
    rectangular_quadmesh.h
    quarter_circle_sector_mesh.h
    fish_mesh.h
    one_d_lagrangian_mesh.h
    two_layer_spine_mesh.h
    single_layer_spine_mesh.h
    bretherton_spine_mesh.h
    quarter_tube_mesh.h
    eighth_sphere_mesh.h
    channel_spine_mesh.h
    collapsible_channel_mesh.h
    triangle_mesh.h
    geompack_mesh.h
    tetgen_mesh.h
    simple_rectangular_tri_mesh.h
    simple_cubic_tet_mesh.h
    simple_cubic_mesh.h
    single_layer_cubic_spine_mesh.h
    channel_with_leaflet_mesh.h
    cylinder_with_flag_mesh.h
    fsi_driven_cavity_mesh.h
    hermite_element_quad_mesh.h
    circular_shell_mesh.h
    tube_mesh.h
    backward_step_mesh.h
    full_circle_mesh.h
    horizontal_single_layer_spine_mesh.h
    xda_tet_mesh.h
    thin_layer_brick_on_tet_mesh.h
    brick_from_tet_mesh.h
    quarter_pipe_mesh.h
    annular_mesh.h
    rectangle_with_hole_mesh.h
    quad_from_triangle_mesh.h
    refineable_tetgen_mesh.h
    gmsh_tet_mesh.h
    rectangle_with_moving_cylinder_mesh.h
    extruded_cube_mesh_from_quad_mesh_with_macro_elements.h)

# Define the template implementation headers (these are not source files as
# we're not trying to build them)
set(TEMPLATE_IMPL_HEADERS
    one_d_mesh.template.cc
    simple_rectangular_quadmesh.template.cc
    rectangular_quadmesh.template.cc
    quarter_circle_sector_mesh.template.cc
    fish_mesh.template.cc
    one_d_lagrangian_mesh.template.cc
    two_layer_spine_mesh.template.cc
    single_layer_spine_mesh.template.cc
    bretherton_spine_mesh.template.cc
    quarter_tube_mesh.template.cc
    eighth_sphere_mesh.template.cc
    channel_spine_mesh.template.cc
    collapsible_channel_mesh.template.cc
    triangle_mesh.template.cc
    geompack_mesh.template.cc
    tetgen_mesh.template.cc
    simple_rectangular_tri_mesh.template.cc
    simple_cubic_tet_mesh.template.cc
    simple_cubic_mesh.template.cc
    single_layer_cubic_spine_mesh.template.cc
    channel_with_leaflet_mesh.template.cc
    cylinder_with_flag_mesh.template.cc
    fsi_driven_cavity_mesh.template.cc
    hermite_element_quad_mesh.template.cc
    circular_shell_mesh.template.cc
    tube_mesh.template.cc
    backward_step_mesh.template.cc
    full_circle_mesh.template.cc
    horizontal_single_layer_spine_mesh.template.cc
    xda_tet_mesh.template.cc
    thin_layer_brick_on_tet_mesh.template.cc
    brick_from_tet_mesh.template.cc
    quarter_pipe_mesh.template.cc
    annular_mesh.template.cc
    rectangle_with_hole_mesh.template.cc
    quad_from_triangle_mesh.template.cc
    refineable_tetgen_mesh.template.cc
    gmsh_tet_mesh.template.cc
    rectangle_with_moving_cylinder_mesh.template.cc
    extruded_cube_mesh_from_quad_mesh_with_macro_elements.template.cc)

# Define the non-mesh headers
set(NON_MESH_HEADERS
    collapsible_channel_domain.h
    fish_domain.h
    channel_with_leaflet_domain.h
    quarter_tube_domain.h
    eighth_sphere_domain.h
    quarter_circle_sector_domain.h
    cylinder_with_flag_domain.h
    topologically_rectangular_domain.h
    tube_domain.h
    full_circle_domain.h
    quarter_pipe_domain.h
    annular_domain.h
    rectangle_with_hole_domain.h)

# Set the name of the library
set(LIBNAME meshes)

# Import the CreateGenericLibrary module to handle the library creation
include(OomphLibraryConfig)
oomph_library_config(
  LIBNAME ${LIBNAME}
  HEADERS ${TEMPLATE_HEADERS} ${NON_MESH_HEADERS}
  HEADERS_NO_COMBINE ${TEMPLATE_IMPL_HEADERS}
  LINKLIBS ${LINKLIBS})

# The meshes library relies on the generic library.
target_link_libraries(meshes INTERFACE oomph::generic)

message(VERBOSE "Leaving meshes subdirectory")
# ------------------------------------------------------------------------------
