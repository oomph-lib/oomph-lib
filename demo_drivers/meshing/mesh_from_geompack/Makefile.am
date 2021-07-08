#Include commands common to every Makefile.am that includes self tests
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Executables with self test
check_PROGRAMS= \
mesh_from_geompack_poisson \
convert_geom_file


# THE EXECUTABLE:
#----------------
# Sources the second executable depends on:
mesh_from_geompack_poisson_SOURCES = mesh_from_geompack_poisson.cc

# Note: The following only works if the libraries have been installed! 

# Required libraries: Only the "generic" and "navier_stokes" libraries, 
# which are accessible via the general library directory which 
# we specify with -L. $(FLIBS) get included just in case
# we decide to use a solver that involves fortran sources. 
mesh_from_geompack_poisson_LDADD = -L@libdir@ -lpoisson -lgeneric  $(EXTERNAL_LIBS) $(FLIBS)




# THE EXECUTABLE:
#----------------
# Sources the second executable depends on:
convert_geom_file_SOURCES = convert_geom_file.cc


EXTRA_DIST += \
box_hole.cs2  box_hole.m2  box_hole.mh2  box_hole.rg2 \
quart_circ.cs2  quart_circ.m2  quart_circ.mh2  quart_circ.rg2 \
rect.cs2  rect.m2  rect.mh2  rect.rg2



