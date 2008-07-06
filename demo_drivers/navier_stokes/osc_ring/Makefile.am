#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables
check_PROGRAMS=osc_ring_alg  osc_ring_macro

# Sources for first executable
osc_ring_alg_SOURCES = osc_ring_alg.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
osc_ring_alg_LDADD = -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


# Sources for second executable
osc_ring_macro_SOURCES = osc_ring_macro.cc

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
osc_ring_macro_LDADD = -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
