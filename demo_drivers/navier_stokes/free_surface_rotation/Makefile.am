#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=free_surface_rotation 

# Sources for executable
free_surface_rotation_SOURCES = free_surface_rotation.cc

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
free_surface_rotation_LDADD = -L@libdir@ -lnavier_stokes -lfluid_interface \
                           -lgeneric \
                           $(EXTERNAL_LIBS) $(FLIBS)
