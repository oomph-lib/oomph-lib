#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=plate

#-------------------------------------------------------------------

# Sources for executable
plate_SOURCES = plate.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
plate_LDADD = -L@libdir@ -lshell -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

# -------------------------------------------------------------------
