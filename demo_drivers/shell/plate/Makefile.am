#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=plate unstructured_clamped_square_plate

#-------------------------------------------------------------------

# Sources for executable
plate_SOURCES = plate.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
plate_LDADD = -L@libdir@ -lshell -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#-------------------------------------------------------------------

# Sources for executable
unstructured_clamped_square_plate_SOURCES = unstructured_clamped_square_plate.cc

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
unstructured_clamped_square_plate_LDADD = -L@libdir@ -lshell -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

# -------------------------------------------------------------------
EXTRA_DIST +=  UnitPlate.1.poly UnitPlate.1.ele UnitPlate.1.node


