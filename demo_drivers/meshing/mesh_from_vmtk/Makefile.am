#Include commands common to every Makefile.am that includes self tests
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Executables with self test
check_PROGRAMS= \
create_fluid_and_solid_surface_mesh_from_fluid_xda_mesh


# Sources the second executable depends on:
create_fluid_and_solid_surface_mesh_from_fluid_xda_mesh_SOURCES = \
create_fluid_and_solid_surface_mesh_from_fluid_xda_mesh.cc

# Required libraries
create_fluid_and_solid_surface_mesh_from_fluid_xda_mesh_LDADD =\
 -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


install: create_fluid_and_solid_surface_mesh_from_fluid_xda_mesh
	if [ ! -e ../../../bin/create_fluid_and_solid_surface_mesh_from_fluid_xda_mesh ]; then $(LN_S) `pwd`/create_fluid_and_solid_surface_mesh_from_fluid_xda_mesh ../../../bin ; fi


EXTRA_DIST +=  iliac.xda replies.txt

