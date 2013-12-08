#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=unstructured_clamped_circular_plate 

#------------------------------------------------------------

# Sources for executable
unstructured_clamped_circular_plate_SOURCES = unstructured_clamped_circular_plate.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
unstructured_clamped_circular_plate_LDADD = -L@libdir@ -lshell -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#--------------------------------------------------------------------
EXTRA_DIST += CircularPlate.1.ele CircularPlate.1.node CircularPlate.1.poly


