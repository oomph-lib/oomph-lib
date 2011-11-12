#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=broken_this_demo customising_output

#----------------------------------------------------------------------

# Sources for executable
broken_this_demo_SOURCES = broken_this_demo.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
# NOT NEEDED AND BREAKS STUFF ON JEAN-LUC'S MACHINE
#broken_this_demo_LDADD = -L@libdir@  $(EXTERNAL_LIBS) $(FLIBS)

#----------------------------------------------------------------------

# Sources for executable
customising_output_SOURCES = customising_output.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
customising_output_LDADD = -L@libdir@  -lpoisson -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
