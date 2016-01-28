#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=axisym_vibrating_shell_non_newtonian

# Sources for executable
axisym_vibrating_shell_non_newtonian_SOURCES = axisym_vibrating_shell_non_newtonian.cc overloaded_element_body.h

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
axisym_vibrating_shell_non_newtonian_LDADD = -L@libdir@ -lgeneralised_newtonian_axisym_navier_stokes -laxisym_navier_stokes -lsolid -lconstitutive -lfluid_interface -lmeshes -lrigid_body -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
