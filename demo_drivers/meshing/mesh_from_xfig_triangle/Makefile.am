#Include commands common to every Makefile.am that includes self tests
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Executables with self test
check_PROGRAMS= \
mesh_from_xfig_poisson \
mesh_from_xfig_navier_stokes


# Helper code
noinst_PROGRAMS=fig2poly


# THE EXECUTABLE:
#----------------
# Sources the second executable depends on:
fig2poly_SOURCES = fig2poly.cc


# THE EXECUTABLE:
#----------------
# Sources the second executable depends on:
mesh_from_xfig_poisson_SOURCES = mesh_from_xfig_poisson.cc

# Note: The following only works if the libraries have been installed! 

# Required libraries: Only the "generic" and "navier_stokes" libraries, 
# which are accessible via the general library directory which 
# we specify with -L. $(FLIBS) get included just in case
# we decide to use a solver that involves fortran sources. 
mesh_from_xfig_poisson_LDADD = -L@libdir@  -lpoisson -lgeneric  \
$(EXTERNAL_LIBS) $(FLIBS)

install: fig2poly
	if [ ! -e ../../../bin/fig2poly ]; then $(LN_S) `pwd`/fig2poly ../../../bin ; fi


# THE EXECUTABLE:
#----------------
# Sources the second executable depends on:
mesh_from_xfig_navier_stokes_SOURCES = mesh_from_xfig_navier_stokes.cc

# Note: The following only works if the libraries have been installed! 

# Required libraries: Only the "generic" and "navier_stokes" libraries, 
# which are accessible via the general library directory which 
# we specify with -L. $(FLIBS) get included just in case
# we decide to use a solver that involves fortran sources. 
mesh_from_xfig_navier_stokes_LDADD = -L@libdir@ -lnavier_stokes -lgeneric \
 $(EXTERNAL_LIBS) $(FLIBS)


EXTRA_DIST += \
 oomph_mesh.fig  oomph_mesh.fig.1.ele  oomph_mesh.fig.1.node \
 oomph_mesh.fig.1.poly  oomph_mesh.fig.poly \
 hole.fig  hole.fig.1.ele  hole.fig.1.node  hole.fig.1.poly  hole.fig.poly

