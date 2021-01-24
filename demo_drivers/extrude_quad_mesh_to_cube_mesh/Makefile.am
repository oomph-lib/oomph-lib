#---------------------------------------------------------------------------
# Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers
#---------------------------------------------------------------------------

# Extra files for distribution; find this using 'ls *py *lay *pvsm *bash *map'
EXTRA_DIST += triangle_meshes

#---------------------------------------------------------------------------

# Name of executable
check_PROGRAMS = \
extrude_triangle_generated_mesh \
extrude_with_macro_element_representation

#---------------------------------------------------------------------------

# Sources for first executable
extrude_triangle_generated_mesh_SOURCES = \
extrude_triangle_generated_mesh.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
extrude_triangle_generated_mesh_LDADD = -L@libdir@ \
-lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#---------------------------------------------------------------------------

# Sources for the executable
extrude_with_macro_element_representation_SOURCES = \
extrude_with_macro_element_representation.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
extrude_with_macro_element_representation_LDADD = -L@libdir@ \
-lpoisson -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#---------------------------------------------------------------------------
