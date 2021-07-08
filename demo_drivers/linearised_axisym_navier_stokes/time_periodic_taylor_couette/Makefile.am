#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=time_periodic_taylor_couette

# Sources for executable
time_periodic_taylor_couette_SOURCES = ../self_starting_BDF2_timestepper.h \
                                       time_periodic_taylor_couette.cc validate.sh

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
time_periodic_taylor_couette_LDADD = -L@libdir@ \
                                     -laxisym_navier_stokes \
                                     -llinearised_axisym_navier_stokes \
                                     -lnavier_stokes -lgeneric \
				     $(EXTERNAL_LIBS) $(FLIBS)

# Include path for library headers: All library headers live in 
# the include directory which we specify with -I
# Automake will replace the variable @includedir@ with the actual
# include directory.
AM_CPPFLAGS += -I@includedir@


