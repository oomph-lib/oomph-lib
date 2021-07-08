#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

#Name of executable
check_PROGRAMS=tensioned_string

#Sources for the executable
tensioned_string_SOURCES = tensioned_string.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
tensioned_string_LDADD = -L@libdir@ -lbeam -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
