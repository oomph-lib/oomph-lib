#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables that should run during the self-test
check_PROGRAMS=two_d_womersley

# Sources for executable
two_d_womersley_SOURCES = two_d_womersley.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
two_d_womersley_LDADD = -L@libdir@ -lwomersley -lnavier_stokes -lgeneric \
				$(EXTERNAL_LIBS) $(FLIBS)

