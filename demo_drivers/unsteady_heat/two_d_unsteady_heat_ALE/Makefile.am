#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables that should run during the self-test
check_PROGRAMS= two_d_unsteady_heat_ALE 

# Sources for executable
two_d_unsteady_heat_ALE_SOURCES = two_d_unsteady_heat_ALE.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
two_d_unsteady_heat_ALE_LDADD = -L@libdir@ -lunsteady_heat -lgeneric \
                                     $(EXTERNAL_LIBS) $(FLIBS)


