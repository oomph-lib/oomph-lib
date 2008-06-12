#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=three_d_entry_flow 

# Sources for executable
three_d_entry_flow_SOURCES = three_d_entry_flow.cc

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
three_d_entry_flow_LDADD = -L@libdir@ -lnavier_stokes -lgeneric \
                           $(EXTERNAL_LIBS) $(FLIBS)
