#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=one_d_advection two_d_advection

# Sources for executable
one_d_advection_SOURCES = one_d_advection.cc validate.sh

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
one_d_advection_LDADD = -L@libdir@ -lflux_transport -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


# Sources for executable
two_d_advection_SOURCES = two_d_advection.cc validate.sh

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
two_d_advection_LDADD = -L@libdir@ -lflux_transport -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

