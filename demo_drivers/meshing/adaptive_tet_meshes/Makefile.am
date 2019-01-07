#Include commands common to every Makefile.am that includes self tests
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Executables with self test
check_PROGRAMS= \
geom_obj_with_boundary


if OOMPH_HAS_CGAL 

# Currently tet mesh adaptation only works with cgal based nearest
# sample point finder in solution transfer

check_PROGRAMS+=\
planar_facet_bounded_mesh_from_gmsh  \
planar_facet_bounded_mesh_from_tetgen \
curved_facet_bounded_mesh_from_gmsh  \
curved_facet_bounded_mesh_from_tetgen

endif


# THE EXECUTABLE:
#----------------

# Local sources that code depends on:
geom_obj_with_boundary_SOURCES = geom_obj_with_boundary.cc 

# Required libraries: 
geom_obj_with_boundary_LDADD = -L@libdir@  -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


# THE EXECUTABLE:
#----------------

# Local sources that code depends on:
planar_facet_bounded_mesh_from_gmsh_SOURCES = planar_facet_bounded_tetmesh.cc \
tetmesh_faceted_surfaces.h

# Required libraries: 
planar_facet_bounded_mesh_from_gmsh_LDADD = -L@libdir@  -lpoisson -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


# THE EXECUTABLE:
#----------------

# Local sources that code depends on:
planar_facet_bounded_mesh_from_tetgen_SOURCES = \
planar_facet_bounded_tetmesh.cc \
tetmesh_faceted_surfaces.h 

# Required libraries: 
planar_facet_bounded_mesh_from_tetgen_LDADD = -L@libdir@  -lpoisson -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

# Macro to switch to tetgen
planar_facet_bounded_mesh_from_tetgen_CXXFLAGS=-DDO_TETGEN


# THE EXECUTABLE:
#----------------

# Local sources that code depends on:
curved_facet_bounded_mesh_from_gmsh_SOURCES = curved_facet_bounded_tetmesh.cc \
tetmesh_faceted_surfaces.h 

# Required libraries: 
curved_facet_bounded_mesh_from_gmsh_LDADD = -L@libdir@  -lpoisson -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


# THE EXECUTABLE:
#----------------

# Local sources that code depends on:
curved_facet_bounded_mesh_from_tetgen_SOURCES = \
curved_facet_bounded_tetmesh.cc \
tetmesh_faceted_surfaces.h 

# Required libraries: 
curved_facet_bounded_mesh_from_tetgen_LDADD = -L@libdir@  -lpoisson -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

# Macro to switch to tetgen
curved_facet_bounded_mesh_from_tetgen_CXXFLAGS=-DDO_TETGEN


EXTRA_DIST += warped_disk_with_triad.lay



