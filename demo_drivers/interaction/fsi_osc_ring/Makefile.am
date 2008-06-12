#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables
check_PROGRAMS=osc_ring_alg osc_ring_macro fsi_osc_ring

# Sources for first executable
fsi_osc_ring_SOURCES = fsi_osc_ring.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
fsi_osc_ring_LDADD = -L@libdir@ -lbeam -lnavier_stokes -lgeneric \
                     $(EXTERNAL_LIBS) $(FLIBS)

# Sources for second executable
osc_ring_alg_SOURCES = osc_ring_alg.cc osc_ring_sarah_asymptotics.h

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
osc_ring_alg_LDADD = -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


# Sources for third executable
osc_ring_macro_SOURCES = osc_ring_macro.cc osc_ring_sarah_asymptotics.h

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
osc_ring_macro_LDADD = -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
