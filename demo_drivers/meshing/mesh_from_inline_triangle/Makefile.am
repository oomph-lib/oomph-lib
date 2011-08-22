#Include commands common to every Makefile.am that includes self tests
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Executables with self test
check_PROGRAMS= \
mesh_from_inline_triangle



# THE EXECUTABLE:
#----------------
# Sources the executable depends on:
mesh_from_inline_triangle_SOURCES = mesh_from_inline_triangle.cc

# Note: The following only works if the libraries have been installed! 

# Required libraries: Only the "generic" and "navier_stokes" libraries, 
# which are accessible via the general library directory which 
# we specify with -L. $(FLIBS) get included just in case
# we decide to use a solver that involves fortran sources. 
mesh_from_inline_triangle_LDADD = -L@libdir@ -lrigid_body -lpoisson \
 -lgeneric  $(EXTERNAL_LIBS) $(FLIBS)
