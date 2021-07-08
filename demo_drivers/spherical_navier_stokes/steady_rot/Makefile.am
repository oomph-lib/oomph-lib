#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=steady_rot 

# Sources for executable
steady_rot_SOURCES = steady_rot.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
steady_rot_LDADD = -L@libdir@ -lspherical_navier_stokes -lgeneric \
                      $(EXTERNAL_LIBS) $(FLIBS)

