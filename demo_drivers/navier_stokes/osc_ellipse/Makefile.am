#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables that should run during the self-test
check_PROGRAMS= osc_quarter_ellipse

# Sources for executable
osc_quarter_ellipse_SOURCES = osc_quarter_ellipse.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
osc_quarter_ellipse_LDADD = -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

