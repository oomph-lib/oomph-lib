#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=orr_sommerfeld

# Sources for executable
orr_sommerfeld_SOURCES = orr_sommerfeld.cc validate.sh

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
orr_sommerfeld_LDADD = -L@libdir@ -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

