
#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables that should run during the self-test
check_PROGRAMS=fsi_jacobian_approximation

# Sources for executable
fsi_jacobian_approximation_SOURCES = fsi_jacobian_approximation.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
fsi_jacobian_approximation_LDADD = -L@libdir@ -lbeam -lnavier_stokes \
                                   -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

