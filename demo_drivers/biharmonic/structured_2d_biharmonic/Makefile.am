#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=structured_2d_biharmonic_bellelement 

#------------------------------------------------------------

# Sources for executable
structured_2d_biharmonic_bellelement_SOURCES = structured_2d_biharmonic_bellelement.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
structured_2d_biharmonic_bellelement_LDADD = -L@libdir@ -lshell -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

