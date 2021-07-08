#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=three_d_entry_flow full_tube 

# Sources for executable
three_d_entry_flow_SOURCES = three_d_entry_flow.cc

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
three_d_entry_flow_LDADD = -L@libdir@ -lnavier_stokes -lgeneric \
                           $(EXTERNAL_LIBS) $(FLIBS)
# Sources for executable
full_tube_SOURCES = full_tube.cc

# Required libraries: 
# $(FLIBS) is included in case the solver involves fortran sources.
full_tube_LDADD = -L@libdir@ -lnavier_stokes -lgeneric \
                           $(EXTERNAL_LIBS) $(FLIBS)
