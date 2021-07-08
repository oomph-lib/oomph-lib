#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=counter_rotating_disks counter_rotating_disks_ref 

# Sources for executable
counter_rotating_disks_SOURCES = counter_rotating_disks.cc \
  linearised_axisym_navier_stokes_elements.h \
  multi_domain_linearised_axisym_navier_stokes_elements.h 

#  linearised_axisym_navier_stokes_elements.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
counter_rotating_disks_LDADD = -L@libdir@ -laxisym_navier_stokes \
			       -lnavier_stokes -lgeneric \
                               $(EXTERNAL_LIBS) $(FLIBS)

# Sources for executable
counter_rotating_disks_ref_SOURCES = counter_rotating_disks_ref.cc \
  refineable_linearised_axisym_navier_stokes_elements.h \
  linearised_axisym_navier_stokes_elements.h \
  multi_domain_linearised_axisym_navier_stokes_elements.h 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
counter_rotating_disks_ref_LDADD = -L@libdir@ -laxisym_navier_stokes \
			           -lnavier_stokes -lgeneric \
                                   $(EXTERNAL_LIBS) $(FLIBS)


EXTRA_DIST += \
  refineable_linearised_axisym_navier_stokes_elements.cc  \
  linearised_axisym_navier_stokes_elements.cc 
