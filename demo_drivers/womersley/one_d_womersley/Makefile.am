#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables that should run during the self-test
check_PROGRAMS=one_d_womersley 

# Sources for executable
one_d_womersley_SOURCES = one_d_womersley.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
one_d_womersley_LDADD = -L@libdir@ -lwomersley -lnavier_stokes -lgeneric \
				$(EXTERNAL_LIBS) $(FLIBS)

