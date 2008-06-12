#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=airy_cantilever

# Sources for executable
airy_cantilever_SOURCES = airy_cantilever.cc

# Required libraries:
# $(FLIBS) is included in case the  solver involves fortran sources.
airy_cantilever_LDADD = -L@libdir@ -lsolid -lconstitutive \
                        -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
