#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables
check_PROGRAMS=lin_unsteady_ring unsteady_ring

# Sources for first executable
lin_unsteady_ring_SOURCES = lin_unsteady_ring.cc

# Required libraries:
# $(FLIBS) get included in case the solver involves fortran sources.
lin_unsteady_ring_LDADD = -L@libdir@ -lbeam -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

# Sources for second executable:
unsteady_ring_SOURCES = unsteady_ring.cc

# Required libraries: 
# $(FLIBS) get included in case the solver involves fortran sources.
unsteady_ring_LDADD = -L@libdir@ -lbeam -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
