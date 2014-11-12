#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=imr_ode

# Sources for executable
imr_ode_SOURCES = imr_ode.cc validate.sh

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
imr_ode_LDADD = -L@libdir@ -lode -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


