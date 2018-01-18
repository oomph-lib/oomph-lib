#-----------------------------------------------------------------
#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers
#-----------------------------------------------------------------

# Names of executable
check_PROGRAMS=\
vorticity_smoother_validation

# Files that aren't mentioned as sources above but are still supposed
# to be distributed
#EXTRA_DIST = run.bash run_vorticity_convergence_test.bash

#-----------------------------------------------------------------

# Sources for executable
vorticity_smoother_validation_SOURCES = vorticity_smoother_validation.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
vorticity_smoother_validation_LDADD = -L@libdir@ -lnavier_stokes -lgeneric \
$(EXTERNAL_LIBS) $(FLIBS)

