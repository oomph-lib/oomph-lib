#Include commands common to every Makefile.am that includes self tests
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Executables with self test
check_PROGRAMS= \
mesh_from_tetgen_poisson \
mesh_from_tetgen_navier_stokes

# THE EXECUTABLE:
#----------------
# Sources the second executable depends on:
mesh_from_tetgen_navier_stokes_SOURCES = mesh_from_tetgen_navier_stokes.cc 

# Note: The following only works if the libraries have been installed! 

# Required libraries: Only the "generic" and "navier_stokes" libraries, 
# which are accessible via the general library directory which 
# we specify with -L. $(FLIBS) get included just in case
# we decide to use a solver that involves fortran sources. 
mesh_from_tetgen_navier_stokes_LDADD = -L@libdir@ -lnavier_stokes \
-lgeneric  $(EXTERNAL_LIBS) $(FLIBS)



# THE EXECUTABLE:
#----------------
# Sources the second executable depends on:
mesh_from_tetgen_poisson_SOURCES = mesh_from_tetgen_poisson.cc 

# Note: The following only works if the libraries have been installed! 

# Required libraries: Only the "generic" and "navier_stokes" libraries, 
# which are accessible via the general library directory which 
# we specify with -L. $(FLIBS) get included just in case
# we decide to use a solver that involves fortran sources. 
mesh_from_tetgen_poisson_LDADD = -L@libdir@ -lpoisson \
-lgeneric  $(EXTERNAL_LIBS) $(FLIBS)


EXTRA_DIST += cube_hole.1.node cube_hole.1.ele cube_hole.poly cube_hole.1.face


