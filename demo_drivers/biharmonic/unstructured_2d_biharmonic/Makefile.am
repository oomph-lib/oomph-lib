#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=unstructured_2d_biharmonic_bellelement unstructured_2d_biharmonic_curvedelement 

#------------------------------------------------------------

# Sources for executable
unstructured_2d_biharmonic_bellelement_SOURCES = unstructured_2d_biharmonic_bellelement.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
unstructured_2d_biharmonic_bellelement_LDADD = -L@libdir@ -lshell -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#-------------------------------------------------------------

# Sources for executable
unstructured_2d_biharmonic_curvedelement_SOURCES = unstructured_2d_biharmonic_curvedelement.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
unstructured_2d_biharmonic_curvedelement_LDADD = -L@libdir@ -lshell -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#--------------------------------------------
EXTRA_DIST += Circle1.1.ele Circle1.1.node Circle1.1.poly

