#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables that should run during the self-test
check_PROGRAMS=flux_control

# Sources for executable
flux_control_SOURCES = flux_control.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
flux_control_LDADD = -L@libdir@ -lwomersley -lnavier_stokes -lgeneric \
				$(EXTERNAL_LIBS) $(FLIBS)

