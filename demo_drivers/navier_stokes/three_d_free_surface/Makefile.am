#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=3d_free_surface

# Sources for executable
3d_free_surface_SOURCES = 3d_free_surface.cc

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
3d_free_surface_LDADD = -L@libdir@ -lnavier_stokes -lsolid -lconstitutive -lfluid_interface -lgeneric \
                     $(EXTERNAL_LIBS) $(FLIBS)
