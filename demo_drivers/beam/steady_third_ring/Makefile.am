#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=steady_third_ring

# Sources for executable
steady_third_ring_SOURCES = steady_third_ring.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
steady_third_ring_LDADD = -L@libdir@ -lbeam -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
