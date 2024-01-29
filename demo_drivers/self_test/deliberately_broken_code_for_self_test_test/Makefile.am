#Include commands common to every Makefile.am that includes self tests
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=code_with_warning


#---------------------------------------------------------------

# Sources for executable
code_with_warning_SOURCES = code_with_warning.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
code_with_warning_LDADD = -L@libdir@ $(EXTERNAL_LIBS) $(FLIBS)


#---------------------------------------------------------------

EXTRA_DIST+=deliberately_broken_code.cc
