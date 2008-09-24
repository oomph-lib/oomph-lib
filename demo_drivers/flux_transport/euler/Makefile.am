#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=one_d_euler two_d_euler couette ff_step ramp

# Sources for executable
one_d_euler_SOURCES = one_d_euler.cc validate.sh

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
one_d_euler_LDADD = -L@libdir@ -lflux_transport -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


# Sources for executable
two_d_euler_SOURCES = two_d_euler.cc validate.sh

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
two_d_euler_LDADD = -L@libdir@ -lflux_transport -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

# Sources for executable
couette_SOURCES = couette.cc validate.sh

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
couette_LDADD = -L@libdir@ -lflux_transport -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

# Sources for executable
ff_step_SOURCES = ff_step.cc validate.sh

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
ff_step_LDADD = -L@libdir@ -lflux_transport -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

# Sources for executable
ramp_SOURCES = ramp.cc validate.sh

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
ramp_LDADD = -L@libdir@ -lflux_transport -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
