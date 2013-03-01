#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables that should run during the self-test
check_PROGRAMS=two_d_unsteady_heat_melt two_d_unsteady_heat_melt_mixed two_d_unsteady_heat_melt_mixed_unstr

# Sources for executable
two_d_unsteady_heat_melt_SOURCES = two_d_unsteady_heat_melt.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
two_d_unsteady_heat_melt_LDADD = -L@libdir@ -lunsteady_heat -lgeneric \
				$(EXTERNAL_LIBS) $(FLIBS)


# Sources for executable
two_d_unsteady_heat_melt_mixed_SOURCES = two_d_unsteady_heat_melt_mixed.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
two_d_unsteady_heat_melt_mixed_LDADD = -L@libdir@ -lunsteady_heat -lgeneric \
				$(EXTERNAL_LIBS) $(FLIBS)


# Sources for executable
two_d_unsteady_heat_melt_mixed_unstr_SOURCES = two_d_unsteady_heat_melt_mixed_unstr.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
two_d_unsteady_heat_melt_mixed_unstr_LDADD = -L@libdir@ -lunsteady_heat -lgeneric \
				$(EXTERNAL_LIBS) $(FLIBS)


EXTRA_DIST += melt.pvsm plot_it.bash

