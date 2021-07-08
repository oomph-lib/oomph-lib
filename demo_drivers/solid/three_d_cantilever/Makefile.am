#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=three_d_cantilever three_d_cantilever_adapt 


#-----------------------------------------------------------------------

# Sources for executable
three_d_cantilever_adapt_SOURCES = three_d_cantilever.cc 

# Required libraries:
# $(FLIBS) is included in case the  solver involves fortran sources.
three_d_cantilever_adapt_LDADD = -L@libdir@ -lsolid -lconstitutive \
                        -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

three_d_cantilever_adapt_CXXFLAGS=-DREFINE


#-----------------------------------------------------------------------

# Sources for executable
three_d_cantilever_SOURCES = three_d_cantilever.cc 

# Required libraries:
# $(FLIBS) is included in case the  solver involves fortran sources.
three_d_cantilever_LDADD = -L@libdir@ -lsolid -lconstitutive \
                        -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#------------------------------------------------------------------------

EXTRA_DIST += comp.mcr