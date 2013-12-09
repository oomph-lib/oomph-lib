#Include commands common to every Makefile.am that includes self tests
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Executables with self test
check_PROGRAMS= \
mesh_from_triangle_poisson \
mesh_from_triangle_poisson_adapt \
mesh_from_triangle_navier_stokes \
mesh_from_triangle_navier_stokes_adapt 



# THE EXECUTABLE:
#----------------
# Sources the second executable depends on:
mesh_from_triangle_poisson_SOURCES = mesh_from_triangle_poisson.cc

# Note: The following only works if the libraries have been installed! 

# Required libraries: Only the "generic" and "navier_stokes" libraries, 
# which are accessible via the general library directory which 
# we specify with -L. $(FLIBS) get included just in case
# we decide to use a solver that involves fortran sources. 
mesh_from_triangle_poisson_LDADD = -L@libdir@ -lpoisson \
 -lgeneric  $(EXTERNAL_LIBS) $(FLIBS)

# THE EXECUTABLE:
#----------------
# Sources the second executable depends on:
mesh_from_triangle_poisson_adapt_SOURCES = mesh_from_triangle_poisson.cc

# Note: The following only works if the libraries have been installed! 

# Required libraries: Only the "generic" and "navier_stokes" libraries, 
# which are accessible via the general library directory which 
# we specify with -L. $(FLIBS) get included just in case
# we decide to use a solver that involves fortran sources. 
mesh_from_triangle_poisson_adapt_LDADD = -L@libdir@ -lpoisson \
 -lgeneric  $(EXTERNAL_LIBS) $(FLIBS)

mesh_from_triangle_poisson_adapt_CXXFLAGS = -DADAPT


# THE EXECUTABLE:
#----------------
# Sources the second executable depends on:
mesh_from_triangle_navier_stokes_SOURCES = mesh_from_triangle_navier_stokes.cc 

# Note: The following only works if the libraries have been installed! 

# Required libraries: Only the "generic" and "navier_stokes" libraries, 
# which are accessible via the general library directory which 
# we specify with -L. $(FLIBS) get included just in case
# we decide to use a solver that involves fortran sources. 
mesh_from_triangle_navier_stokes_LDADD = -L@libdir@ -lnavier_stokes -lgeneric \
 $(EXTERNAL_LIBS) $(FLIBS)

# THE EXECUTABLE:
#----------------
# Sources the second executable depends on:
mesh_from_triangle_navier_stokes_adapt_SOURCES = mesh_from_triangle_navier_stokes.cc 

# Note: The following only works if the libraries have been installed! 

# Required libraries: Only the "generic" and "navier_stokes" libraries, 
# which are accessible via the general library directory which 
# we specify with -L. $(FLIBS) get included just in case
# we decide to use a solver that involves fortran sources. 
mesh_from_triangle_navier_stokes_adapt_LDADD = -L@libdir@ -lnavier_stokes -lgeneric \
 $(EXTERNAL_LIBS) $(FLIBS)

mesh_from_triangle_navier_stokes_adapt_CXXFLAGS = -DADAPT

EXTRA_DIST += box_hole.poly box_hole.1.node box_hole.1.ele box_hole.1.poly \
flow_past_box.poly flow_past_box.1.node flow_past_box.1.ele \
flow_past_box.1.poly box_hole.1.edge flow_past_box.1.edge \
box_hole_adapt.poly box_hole_adapt.1.node box_hole_adapt.1.ele \
box_hole_adapt.1.poly flow_past_box_adapt.poly flow_past_box_adapt.1.node \
flow_past_box_adapt.1.ele flow_past_box_adapt.1.poly 

