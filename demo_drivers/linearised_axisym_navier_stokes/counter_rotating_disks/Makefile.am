#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=counter_rotating_disks counter_rotating_disks_ref

# Sources for executable
counter_rotating_disks_SOURCES = ../self_starting_BDF2_timestepper.h \
                                 counter_rotating_disks.cc validate.sh

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
counter_rotating_disks_LDADD = -L@libdir@ -llinearised_axisym_navier_stokes \
                               -laxisym_navier_stokes \
                               -lnavier_stokes -lgeneric \
                               $(EXTERNAL_LIBS) $(FLIBS)

# Sources for executable
counter_rotating_disks_ref_SOURCES = ../self_starting_BDF2_timestepper.h \
                                     counter_rotating_disks_ref.cc validate.sh

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
counter_rotating_disks_ref_LDADD = -L@libdir@ \
				   -llinearised_axisym_navier_stokes \
                                   -laxisym_navier_stokes \
                                   -lnavier_stokes -lgeneric \
                                   $(EXTERNAL_LIBS) $(FLIBS)

# Include path for library headers: All library headers live in 
# the include directory which we specify with -I
# Automake will replace the variable @includedir@ with the actual
# include directory.
AM_CPPFLAGS += -I@includedir@

