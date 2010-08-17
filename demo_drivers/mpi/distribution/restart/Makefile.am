#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables that should run during the self-test
check_PROGRAMS=two_d_unsteady_heat_2adapt_load_balance

#----------------------------------------------

# Sources for executable
two_d_unsteady_heat_2adapt_load_balance_SOURCES = two_d_unsteady_heat_2adapt_load_balance.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
two_d_unsteady_heat_2adapt_load_balance_LDADD = -L@libdir@ -lunsteady_heat -lgeneric \
				    $(EXTERNAL_LIBS) $(FLIBS)

#----------------------------------------------

EXTRA_DIST+=doubly_adaptive_partitioning_load_balance.dat \
            load_balance.mcr  show_dist.mcr
